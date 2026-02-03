import streamlit as st
import sqlite3
import pandas as pd
import json
import re
import requests
from Bio import Entrez
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from openai import OpenAI
import time
import db_handler as db  # [중요] DB 핸들러 임포트

# ===================== [앱 시작 시 DB 동기화] =====================
# 앱이 켜질 때 GitHub에서 최신 DB 파일을 받아옵니다.
if 'db_synced' not in st.session_state:
    try:
        with st.spinner("데이터 동기화 중..."):
            db.pull_db()
        st.session_state.db_synced = True
    except Exception as e:
        st.warning(f"DB 동기화 실패: {e}")

# ===================== [설정 및 초기화] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터 Pro", layout="wide", page_icon="🏥")

with st.sidebar:
    st.header("⚙️ 기본 설정")
    openai_api_key = st.text_input("OpenAI API Key", type="password")
    email_address = st.text_input("Email (PubMed용)", value="your_email@example.com")
    
    st.divider()
    st.header("📢 텔레그램 설정")
    st.caption("단톡방 ID는 보통 마이너스(-)로 시작합니다.")
    telegram_token = st.text_input("Bot Token", type="password")
    chat_id = st.text_input("Chat ID")

Entrez.email = email_address
DB_NAME = 'kmd_papers_v5_column.db' 

# ===================== [1. DB 관리] =====================
def init_db():
    # db_handler에서 처리하므로 여기선 생략 가능하지만, 안전장치로 둠
    pass

def get_papers_by_date(target_date_str):
    conn = sqlite3.connect(DB_NAME)
    try:
        query = "SELECT * FROM papers WHERE date_published = ?"
        df = pd.read_sql(query, conn, params=(target_date_str,))
    except:
        df = pd.DataFrame()
    conn.close()
    return df

def get_daily_column(date_str):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT content FROM daily_columns WHERE date_id = ?", (date_str,))
        result = cursor.fetchone()
        return result[0] if result else None
    except: return None
    finally: conn.close()

def save_daily_column(date_str, content):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO daily_columns VALUES (?, ?, ?)", 
                   (date_str, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()
    # [중요] 저장 후 GitHub 업로드
    db.push_db()

def get_blog_post(date_str, target_type):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT content FROM blog_posts WHERE date_id = ? AND target_type = ?", (date_str, target_type))
        result = cursor.fetchone()
        return result[0] if result else None
    except: return None
    finally: conn.close()

def save_blog_post(date_str, target_type, content):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO blog_posts VALUES (?, ?, ?, ?)", 
                   (date_str, target_type, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()
    # [중요] 저장 후 GitHub 업로드
    db.push_db()

def delete_papers(pmid_list):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    if pmid_list:
        placeholders = ', '.join('?' for _ in pmid_list)
        cursor.execute(f"DELETE FROM papers WHERE pmid IN ({placeholders})", pmid_list)
        conn.commit()
    conn.close()
    # [중요] 삭제 후 GitHub 업로드
    db.push_db()

def check_if_exists(pmid):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,))
    exists = cursor.fetchone() is not None
    conn.close()
    return exists

# ===================== [2. Full Text Fetcher] =====================
def fetch_pmc_fulltext(pmid):
    try:
        link_results = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        if not link_results or not link_results[0]['LinkSetDb']:
            return None, "Abstract Only"
        pmc_id = link_results[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        xml_data = handle.read()
        root = ET.fromstring(xml_data)
        full_text = "".join([text for body in root.findall(".//body") for text in body.itertext()])
        return full_text[:25000] if len(full_text) > 500 else None, "✅ Full Text (PMC)"
    except Exception as e:
        return None, f"Error: {str(e)}"

# ===================== [3. AI 분석 로직] =====================
def analyze_paper_strict(paper_data, api_key):
    client = OpenAI(api_key=api_key)
    prompt = f"""
    너는 임상 한의학 논문 분류 전문가다.
    
    [필수 규칙 1: 중재법 분류]
    - 침, 뜸, 부항, 한약, 약침, 추나 중 택1. 해당 없으면 "기타".

    [필수 규칙 2: 신체부위 분류]
    - 두경부, 척추/허리, 상지, 하지, 내장기/전신 중 택1.

    [JSON 형식]
    {{
        "korean_title": "한글 제목",
        "study_design": "연구 유형",
        "intervention_category": "카테고리",
        "target_body_part": "신체부위",
        "specific_point": "상세 중재 내용",
        "clinical_score": 8,
        "summary": "3줄 요약",
        "icd_code": "코드",
        "full_text_status": "Abstract Check"
    }}

    Title: {paper_data['title']}
    Abstract: {paper_data['abstract']}
    """
    try:
        response = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0
        )
        data = json.loads(re.search(r'\{.*\}', response.choices[0].message.content.strip(), re.DOTALL).group())
        if "DROP" in str(data.get("study_design", "")): return {"error": "DROP: 임상 연구 아님"}
        return data
    except Exception as e:
        return {"error": str(e)}

# ===================== [4. PubMed 검색] =====================
def simple_keyword_classify(text):
    text = text.lower()
    if "acupuncture" in text or "needling" in text: return "침"
    elif "moxibustion" in text: return "뜸"
    elif "cupping" in text: return "부항"
    elif "herbal" in text or "decoction" in text: return "한약"
    elif "pharmacopuncture" in text or "injection" in text: return "약침"
    elif "chuna" in text or "manipulation" in text: return "추나"
    else: return "기타"

def search_pubmed_raw(start_date, end_date, max_results):
    str_start = start_date.strftime("%Y/%m/%d")
    str_end = end_date.strftime("%Y/%m/%d")
    search_term = """
    ("TCM" OR "Traditional chinese medicine" OR "Herbal medicine" OR "Medicine, Korean Traditional" 
    OR "Acupuncture" OR "Moxibustion" OR "Cupping Therapy" OR "Pharmacopuncture" OR "Chuna") 
    AND (hasabstract[text]) AND ("Humans"[Mesh]) 
    AND ("Case Reports"[ptyp] OR "Clinical Trial"[ptyp] OR "Randomized Controlled Trial"[ptyp] OR "Systematic Review"[ptyp] OR "Cohort Studies"[Mesh])
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, mindate=str_start, maxdate=str_end, datetype="pdat", retmax=max_results)
        record = Entrez.read(handle)
        id_list = record["IdList"]
    except: return []

    if not id_list: return []

    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
    except: return []

    raw_papers = []
    for article in records['PubmedArticle']:
        try:
            pmid = str(article['MedlineCitation']['PMID'])
            title = article['MedlineCitation']['Article']['ArticleTitle']
            abstract_list = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
            abstract = " ".join(abstract_list) if abstract_list else ""
            
            raw_papers.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract,
                "predicted_category": simple_keyword_classify(title + abstract),
                "is_saved": check_if_exists(pmid)
            })
        except: continue
    return raw_papers

# ===================== [5. 데일리 브리핑 생성기] =====================
def generate_daily_briefing_pro_v3(date_str, papers_df, api_key, model_choice):
    client = OpenAI(api_key=api_key)
    top_papers = papers_df.sort_values(by='clinical_score', ascending=False).head(10)
    
    if top_papers.empty: return "분석할 논문이 없습니다."

    analyzed_data = []
    prog_bar = st.progress(0)
    status_text = st.empty()

    for idx, (_, row) in enumerate(top_papers.iterrows()):
        prog_bar.progress((idx+1)/len(top_papers))
        status_text.text(f"🔍 안전 분석 모드 동작 중... ({idx+1}): {row['title_kr']}")
        
        full_text, ft_status = fetch_pmc_fulltext(row['pmid'])
        content_source = full_text if full_text else row['abstract']
        
        pico_prompt = f"""
        이 논문을 PICO 구조로 분석하라.
        [규칙] Full Name 변환 및 수치 정보 포함.
        Title: {row['title_kr']}
        Text: {content_source[:15000]}
        """
        try:
            pico_res_text = client.chat.completions.create(
                model="gpt-4o-mini", 
                messages=[{"role": "user", "content": pico_prompt}],
                temperature=0.0
            ).choices[0].message.content
        except: pico_res_text = "분석 실패"

        analyzed_data.append({
            "pmid": row['pmid'],
            "title": row['title_kr'],
            "score": row['clinical_score'],
            "study_design": row['study_design'],
            "source": ft_status,
            "detail_analysis": pico_res_text
        })

    status_text.empty()
    prog_bar.empty()

    final_prompt = f"""
    당신은 한의학 에디터입니다. 상위 7개(Pick 2 + News 5) 논문 브리핑을 작성하세요.
    [필수] 원문 링크 포함: `🔗 원문: https://pubmed.ncbi.nlm.nih.gov/[PMID]`
    
    [출력 포맷]
    📅 **{date_str} 한의 임상 브리핑**
    
    🥇 **Today's Pick 1: [제목]**
    ([연구유형] / ⭐[점수])
    - 🎯 **Point:** ...
    - 💊 **Method:** ...
    - 📊 **Result:** ...
    - 🔗 **원문:** https://pubmed.ncbi.nlm.nih.gov/[PMID]
    
    [입력 데이터]
    {json.dumps(analyzed_data, ensure_ascii=False)}
    """
    
    try:
        if "o1" in model_choice:
            return client.chat.completions.create(model=model_choice, messages=[{"role": "user", "content": final_prompt}]).choices[0].message.content
        else:
            return client.chat.completions.create(model=model_choice, messages=[{"role": "user", "content": final_prompt}], temperature=0.7).choices[0].message.content
    except Exception as e: return f"생성 실패: {e}"

# ===================== [6. 블로그 아티클 생성기] =====================
def generate_blog_article(date_str, papers_df, api_key, model_choice, target_audience="doctor"):
    client = OpenAI(api_key=api_key)
    top_paper = papers_df.sort_values(by='clinical_score', ascending=False).iloc[0]
    full_text, ft_status = fetch_pmc_fulltext(top_paper['pmid'])
    content_source = full_text if full_text else top_paper['abstract']
    
    prompt = f"""
    당신은 전문 의학 블로거입니다. 이 논문으로 블로그 글을 작성하세요.
    타겟: {'전문가(한의사)' if target_audience == 'doctor' else '일반 환자'}
    
    [논문 정보]
    제목: {top_paper['title_kr']}
    PMID: {top_paper['pmid']}
    내용: {content_source[:20000]}
    """
    
    try:
        if "o1" in model_choice:
            return client.chat.completions.create(model=model_choice, messages=[{"role": "user", "content": prompt}]).choices[0].message.content
        else:
            return client.chat.completions.create(model=model_choice, messages=[{"role": "user", "content": prompt}], temperature=0.7).choices[0].message.content
    except Exception as e: return f"블로그 생성 실패: {e}"

# ===================== [7. UI 구성] =====================
st.title("🏥 한의학 논문 AI 큐레이터 Pro")
st.markdown("---")

tab_briefing, tab_blog, tab_archive, tab_search = st.tabs(["📝 데일리 브리핑", "✍️ 블로그/수익화", "📚 보관함", "🔎 검색"])

# --- [Tab 1: 데일리 브리핑 & 텔레그램] ---
with tab_briefing:
    c1, c2 = st.columns([1, 2])
    with c1:
        st.subheader("🗓 브리핑 설정")
        target_date = st.date_input("날짜", value=datetime.now())
        target_date_str = target_date.strftime("%Y-%m-%d")
        daily_papers = get_papers_by_date(target_date_str)
        st.info(f"논문 수: {len(daily_papers)}건")
        
        model_option = st.radio("AI 모델:", ["gpt-4o", "o1-preview", "gpt-4o-mini"], index=0)
        
        if st.button("✨ 안전 모드 브리핑 생성"):
            if daily_papers.empty: st.error("논문 없음")
            elif not openai_api_key: st.error("Key 없음")
            else:
                briefing = generate_daily_briefing_pro_v3(target_date_str, daily_papers, openai_api_key, model_option)
                save_daily_column(target_date_str, briefing)
                st.success("브리핑 생성 및 저장 완료!") # GitHub 업로드 됨
                st.rerun()

    with c2:
        st.subheader("📨 공유 및 전송")
        content = get_daily_column(target_date_str)
        
        if content:
            st.markdown("##### 🚀 텔레그램 전송")
            user_footer = st.text_area("📢 추가 코멘트", height=70)
            final_msg = content
            if user_footer: final_msg += f"\n\n--------------------------------\n📢 **Editor's Note**\n{user_footer}"

            if st.button("✈️ 텔레그램 전송", type="primary"):
                if not telegram_token or not chat_id:
                    st.error("설정에서 토큰/ID를 입력하세요.")
                else:
                    try:
                        url = f"https://api.telegram.org/bot{telegram_token}/sendMessage"
                        res = requests.post(url, json={"chat_id": chat_id, "text": final_msg, "parse_mode": "Markdown"})
                        if res.status_code == 200: st.toast("전송 성공!"); st.success("전송 완료")
                        else: st.error(f"실패: {res.text}")
                    except Exception as e: st.error(f"에러: {e}")

            st.divider()
            st.code(final_msg, language='markdown')
            with st.expander("미리보기"): st.markdown(final_msg)
        else:
            st.warning("생성된 브리핑이 없습니다.")

# --- [Tab 2: 블로그] ---
with tab_blog:
    c_b1, c_b2 = st.columns([1, 3])
    with c_b1:
        st.subheader("✒️ 블로그 생성")
        b_date = st.date_input("날짜", value=datetime.now(), key="blog_date")
        b_date_str = b_date.strftime("%Y-%m-%d")
        b_papers = get_papers_by_date(b_date_str)
        st.info(f"후보: {len(b_papers)}건")
        b_model = st.selectbox("모델:", ["gpt-4o", "o1-preview"], index=0)
        target_type = st.radio("타겟:", ["👨‍⚕️ 전문가용", "😊 환자용"])
        
        if st.button("✍️ 글 생성"):
            if b_papers.empty: st.error("논문 없음")
            elif not openai_api_key: st.error("Key 없음")
            else:
                with st.spinner("작성 중..."):
                    article = generate_blog_article(b_date_str, b_papers, openai_api_key, b_model, "doctor" if "전문가" in target_type else "patient")
                    save_blog_post(b_date_str, "doctor" if "전문가" in target_type else "patient", article)
                    st.success("완료! (GitHub 자동 저장됨)")
                    st.rerun()

    with c_b2:
        st.subheader("📄 미리보기")
        t1, t2 = st.tabs(["👨‍⚕️ 전문가용", "😊 환자용"])
        with t1:
            post = get_blog_post(b_date_str, "doctor")
            if post: st.markdown(post); st.divider(); st.code(post)
            else: st.info("없음")
        with t2:
            post = get_blog_post(b_date_str, "patient")
            if post: st.markdown(post); st.divider(); st.code(post)
            else: st.info("없음")

# --- [Tab 3: 보관함 (필터 적용)] ---
with tab_archive:
    # 1. DB에서 데이터 가져오기
    df_all = pd.read_sql("SELECT * FROM papers", sqlite3.connect(DB_NAME))
    
    if df_all.empty:
        st.info("보관함이 비어있습니다.")
    else:
        st.subheader("🔍 필터링")
        
        # 2. 필터 UI
        cats = sorted(df_all['intervention_category'].unique().tolist())
        sel_cats = st.multiselect("중재법 선택", cats, default=cats)
        
        if 'archive_body_part' not in st.session_state: st.session_state.archive_body_part = "전체"
        def btn_col(part): return "primary" if st.session_state.archive_body_part == part else "secondary"
        
        parts = ["두경부", "척추/허리", "상지", "하지", "내장기/전신", "전체"]
        cols = st.columns(6)
        for i, part in enumerate(parts):
            if cols[i].button(part, key=f"p_{i}", type=btn_col(part), use_container_width=True):
                st.session_state.archive_body_part = part
                st.rerun()

        # 3. 데이터 필터링
        df_filt = df_all.copy()
        if sel_cats: df_filt = df_filt[df_filt['intervention_category'].isin(sel_cats)]
        if st.session_state.archive_body_part != "전체":
            df_filt = df_filt[df_filt['target_body_part'] == st.session_state.archive_body_part]

        st.divider()
        st.subheader(f"📚 목록 ({len(df_filt)}건)")
        
        if not df_filt.empty:
            # 삭제 체크박스 및 URL 링크 생성
            df_filt.insert(0, "del", False)
            df_filt["url"] = "https://pubmed.ncbi.nlm.nih.gov/" + df_filt["pmid"]
            
            # 4. 데이터 에디터 표시 (여기에 'date_published' 추가됨!)
            edited = st.data_editor(
                df_filt,
                column_config={
                    "del": st.column_config.CheckboxColumn("삭제", width="small"),
                    "url": st.column_config.LinkColumn("Link", display_text="🔗", width="small"),
                    "date_published": st.column_config.TextColumn("수집일", width="small"), # [복구됨]
                    "title_kr": st.column_config.TextColumn("제목", width="large"),
                    "target_body_part": st.column_config.TextColumn("부위", width="small"),
                    "intervention_category": st.column_config.TextColumn("중재", width="small"),
                    "clinical_score": st.column_config.NumberColumn("점수", format="%d점"),
                },
                # 컬럼 순서 지정 (수집일을 앞쪽으로 배치)
                column_order=["del", "url", "date_published", "clinical_score", "intervention_category", "target_body_part", "title_kr", "summary"],
                hide_index=True, 
                use_container_width=True
            )
            
            # 5. 삭제 로직
            if st.button("🗑️ 삭제 확인"):
                to_del = edited[edited["del"]]['pmid'].tolist()
                if to_del:
                    delete_papers(to_del)
                    st.success("삭제됨 (GitHub 자동 동기화)")
                    st.rerun()
        else: st.warning("조건에 맞는 논문이 없습니다.")

# --- [Tab 4: 검색] ---
with tab_search:
    c1, c2 = st.columns(2)
    with c1: s_date = st.date_input("시작", value=datetime.now()-timedelta(days=2))
    with c2: e_date = st.date_input("종료", value=datetime.now())
    limit = st.slider("개수", 10, 100, 50)
    
    if 'search_res' not in st.session_state: st.session_state.search_res = None
    if st.button("1. 검색"):
        with st.spinner(".."): st.session_state.search_res = search_pubmed_raw(s_date, e_date, limit)
        
    if st.session_state.search_res:
        df = pd.DataFrame(st.session_state.search_res)
        df.insert(0, "Sel", ~df['is_saved'])
        edited = st.data_editor(df, column_config={"Sel": st.column_config.CheckboxColumn("선택")}, hide_index=True)
        targets = edited[edited["Sel"]]
        
        if st.button(f"2. {len(targets)}건 분석 및 저장"):
            if not openai_api_key: st.error("Key 없음")
            else:
                conn = sqlite3.connect(DB_NAME); cur = conn.cursor()
                bar = st.progress(0)
                full_list = [p for p in st.session_state.search_res if p['pmid'] in targets['pmid'].tolist()]
                
                # DB 테이블 생성 체크
                cur.execute('''CREATE TABLE IF NOT EXISTS papers (pmid TEXT PRIMARY KEY, date_published TEXT, title_kr TEXT, intervention_category TEXT, target_body_part TEXT, specific_point TEXT, study_design TEXT, clinical_score INTEGER, summary TEXT, original_title TEXT, abstract TEXT, icd_code TEXT, full_text_status TEXT)''')
                
                for i, p in enumerate(full_list):
                    bar.progress((i+1)/len(full_list))
                    res = analyze_paper_strict(p, openai_api_key)
                    if "error" not in res:
                        cur.execute('INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', (
                            p['pmid'], datetime.now().strftime('%Y-%m-%d'),
                            res.get('korean_title'), res.get('intervention_category'),
                            res.get('target_body_part'), res.get('specific_point'),
                            res.get('study_design'), res.get('clinical_score'),
                            res.get('summary'), p['title'], p['abstract'], 
                            res.get('icd_code'), "Abstract Saved"
                        ))
                        conn.commit()
                conn.close()
                # [중요] 저장 완료 후 GitHub 업로드
                db.push_db()
                
                st.success("분석 및 저장 완료! (GitHub 동기화 됨)")
                st.session_state.search_res = None
                time.sleep(1)
                st.rerun()

# 앱 실행 시 DB 체크 (임포트 시 실행되지만 안전상 한 번 더)
if __name__ == "__main__":
    if not st.session_state.get('db_synced'):
        db.pull_db()
        st.session_state.db_synced = True


