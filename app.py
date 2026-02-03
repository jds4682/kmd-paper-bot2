import streamlit as st
import sqlite3
import pandas as pd
import json
import re
from Bio import Entrez
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from openai import OpenAI
import time

# ===================== [설정 및 초기화] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터 Pro", layout="wide", page_icon="🏥")

with st.sidebar:
    st.header("⚙️ 설정")
    openai_api_key = st.text_input("OpenAI API Key", type="password")
    email_address = st.text_input("Email (PubMed용)", value="your_email@example.com")
    st.info("💡 Top 논문은 원문(Full Text)을 확보하여 심층 분석합니다.")

Entrez.email = email_address
DB_NAME = 'kmd_papers_v5_column.db' 

# ===================== [1. DB 관리 (자동 업데이트 기능 추가)] =====================
def init_db():
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    
    # 1. 논문 테이블 생성 (기본)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS papers (
            pmid TEXT PRIMARY KEY,
            date_published TEXT,
            title_kr TEXT,
            intervention_category TEXT, 
            target_body_part TEXT,      
            specific_point TEXT,        
            study_design TEXT,
            clinical_score INTEGER,
            summary TEXT,
            original_title TEXT,
            abstract TEXT,
            icd_code TEXT,
            full_text_status TEXT
        )
    ''')
    
    # 🚨 DB 마이그레이션 (에러 해결 핵심 로직)
    # 기존 DB에 'full_text_status' 컬럼이 없으면 추가해주는 코드
    cursor.execute("PRAGMA table_info(papers)")
    columns = [info[1] for info in cursor.fetchall()]
    if 'full_text_status' not in columns:
        try:
            cursor.execute("ALTER TABLE papers ADD COLUMN full_text_status TEXT")
            st.toast("시스템: DB 업데이트 완료 (full_text_status 컬럼 추가됨)")
        except:
            pass # 이미 있으면 패스

    # 2. 데일리 브리핑 테이블
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS daily_columns (
            date_id TEXT PRIMARY KEY,
            content TEXT,
            created_at TEXT
        )
    ''')
    
    # 3. 블로그 포스트 테이블 (구버전 충돌 방지)
    try:
        cursor.execute("SELECT target_type FROM blog_posts LIMIT 1")
    except sqlite3.OperationalError:
        cursor.execute("DROP TABLE IF EXISTS blog_posts")
        conn.commit()

    cursor.execute('''
        CREATE TABLE IF NOT EXISTS blog_posts (
            date_id TEXT,
            target_type TEXT, 
            content TEXT,
            created_at TEXT,
            PRIMARY KEY (date_id, target_type)
        )
    ''')
    
    # 4. [NEW] 시스템 설정 테이블 (자동화 제어용)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS system_config (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    ''')
    
    conn.commit()
    conn.close()

# 설정 저장/로드 함수
def set_config(key, value):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute("INSERT OR REPLACE INTO system_config VALUES (?, ?)", (key, str(value)))
    conn.commit()
    conn.close()

def get_config(key):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute("SELECT value FROM system_config WHERE key=?", (key,))
    res = cur.fetchone()
    conn.close()
    return res[0] if res else "False"

def get_papers_by_date(target_date_str):
    conn = sqlite3.connect(DB_NAME)
    try:
        query = "SELECT * FROM papers WHERE date_published = ?"
        df = pd.read_sql(query, conn, params=(target_date_str,))
    except:
        df = pd.DataFrame()
    conn.close()
    return df

# 브리핑 저장/로드
def get_daily_column(date_str):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("SELECT content FROM daily_columns WHERE date_id = ?", (date_str,))
    result = cursor.fetchone()
    conn.close()
    return result[0] if result else None

def save_daily_column(date_str, content):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO daily_columns VALUES (?, ?, ?)", 
                   (date_str, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()

# 블로그 저장/로드
def get_blog_post(date_str, target_type):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT content FROM blog_posts WHERE date_id = ? AND target_type = ?", (date_str, target_type))
        result = cursor.fetchone()
        return result[0] if result else None
    except sqlite3.OperationalError:
        return None

def save_blog_post(date_str, target_type, content):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO blog_posts VALUES (?, ?, ?, ?)", 
                   (date_str, target_type, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()

def delete_papers(pmid_list):
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    if pmid_list:
        placeholders = ', '.join('?' for _ in pmid_list)
        cursor.execute(f"DELETE FROM papers WHERE pmid IN ({placeholders})", pmid_list)
        conn.commit()
    conn.close()

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
        handle.close()
        
        root = ET.fromstring(xml_data)
        full_text = ""
        for body in root.findall(".//body"):
            for text in body.itertext():
                full_text += text + " "
        
        if len(full_text) > 500:
            return full_text[:25000], "✅ Full Text (PMC)" 
        else:
            return None, "PMC XML Empty"
    except Exception as e:
        return None, f"Error: {str(e)}"

# ===================== [3. AI 분석 로직] =====================
def analyze_paper_strict(paper_data, api_key):
    client = OpenAI(api_key=api_key)
    prompt = f"""
    너는 임상 한의학 논문 심사관이다.
    
    [필수 규칙]
    0. **동물/세포 실험은 무조건 DROP.**
    1. 'clinical_score': 1~10점 (근골격/소화기/통증 등 로컬 다빈도 질환 가산점).
    2. 'specific_point': 처방명(구성/g수), 혈자리 필수.
    3. 'study_design': RCT, SR, Case Report, Cohort.

    [JSON 형식]
    {{
        "korean_title": "한글 제목",
        "study_design": "연구 유형",
        "intervention_category": "카테고리",
        "specific_point": "상세 중재 내용",
        "target_body_part": "신체부위",
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
    if "acupuncture" in text or "needling" in text: return "1_침구치료"
    elif "herbal" in text or "decoction" in text: return "2_한약치료"
    elif "chuna" in text or "manipulation" in text: return "5_추나/도수"
    else: return "7_기타/복합"

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
        handle.close()
        id_list = record["IdList"]
    except: return []

    if not id_list: return []

    try:
        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        handle.close()
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

        [🚨 안전 분석 규칙]
        1. **약어(Acronym):** 텍스트 내에서 정의된 Full Name을 찾아라. 정의가 없으면 'Unknown'으로 표기.
        2. **상세 정보:** 약재 용량(g), 횟수 등이 없으면 '본문 미기재' 표기.
        
        [분석 항목]
        - valid: true/false
        - P: 환자 정보
        - I: 중재 (약어는 Full name 변환, 없으면 '정보 없음')
        - C: 대조군
        - O: 결과 (p-value 포함, 없으면 '수치 미기재')

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
    당신은 팩트 체크를 중요시하는 한의학 에디터입니다.
    제공된 {len(analyzed_data)}개의 분석 데이터를 검토하여 **상위 7개(Pick 2 + News 5)**만 선별하여 브리핑을 작성하세요.

    [선별 및 작성 규칙]
    1. **엄격한 필터링:** 정보가 부족하거나 약어가 불분명한 논문은 제외하세요.
    2. **내용 충실:** News 섹션도 상세하게 작성.
    3. **필수 링크:** 각 논문의 마지막 줄에 반드시 원문 링크를 달아주세요.
       - 형식: `🔗 원문: https://pubmed.ncbi.nlm.nih.gov/[PMID]`
       - 입력 데이터의 'pmid' 값을 사용하세요.
    
    [출력 포맷]
    📅 **{date_str} 한의 임상 브리핑**
    
    🥇 **Today's Pick 1: [제목]**
    ([연구유형] / ⭐[점수])
    - 🎯 **Point:** ...
    - 💊 **Method:** ...
    - 📊 **Result:** ...
    - 🔎 **Check:** (원문 분석 여부)
    - 🔗 **원문:** https://pubmed.ncbi.nlm.nih.gov/[PMID]

    🥈 **Today's Pick 2: [제목]**
    (동일 양식)

    --------------------------------
    
    📰 **Clinical News Top 5 (상세)**
    
    1️⃣ **[제목]**
       - 📝 **치료:** ...
       - 📉 **결과:** ...
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
    
    if target_audience == "doctor":
        prompt = f"""
        당신은 전문 의학 블로거입니다. 오늘의 Top 논문 1편을 선정하여 **1500자 내외의 심층 전문 아티클**을 작성하세요.

        [목표] '구글 애드센스' 수익형 블로그에 적합한 고품질 콘텐츠.

        [논문 데이터]
        제목: {top_paper['title_kr']} ({top_paper['original_title']})
        PMID: {top_paper['pmid']}
        내용: {content_source[:20000]}
        """
    else:
        prompt = f"""
        당신은 다정하고 실력 있는 동네 한의원 원장님입니다. 
        이 논문의 내용을 바탕으로 일반인 환자들이 읽기 쉬운 **네이버 블로그 포스팅**을 작성하세요.

        [목표] 환자들에게 신뢰를 주고 내원 유도.

        [논문 데이터]
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
init_db()

st.title("🏥 한의학 논문 AI 큐레이터 Pro")
st.markdown("---")

tab_briefing, tab_blog, tab_archive, tab_search = st.tabs(["📝 데일리 브리핑", "✍️ 블로그/수익화", "📚 보관함", "🔎 검색"])

# --- [Tab 1: 데일리 브리핑] ---
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
                st.rerun()

    with c2:
        st.subheader("📨 공유 및 커스텀")
        content = get_daily_column(target_date_str)
        
        if content:
            st.markdown("##### 📢 추가 코멘트")
            user_footer = st.text_area("하단 안내문구", height=100)
            
            final_display_text = content
            if user_footer:
                final_display_text += "\n\n--------------------------------\n📢 **Editor's Note**\n" + user_footer
            
            st.code(final_display_text, language='markdown')
            st.caption("👆 우측 상단 아이콘 클릭 시 전체 복사")
            
            with st.expander("미리보기"):
                st.markdown(final_display_text)
        else:
            st.warning("생성된 브리핑이 없습니다.")

# --- [Tab 2: 블로그/수익화] ---
with tab_blog:
    c_b1, c_b2 = st.columns([1, 3])
    with c_b1:
        st.subheader("✒️ 블로그 생성")
        b_date = st.date_input("블로그 날짜", value=datetime.now(), key="blog_date")
        b_date_str = b_date.strftime("%Y-%m-%d")
        b_papers = get_papers_by_date(b_date_str)
        st.info(f"후보 논문: {len(b_papers)}건")
        
        b_model = st.selectbox("작성 모델:", ["gpt-4o", "o1-preview"], index=0)
        
        st.divider()
        target_type = st.radio(
            "타겟 독자:", 
            ["👨‍⚕️ 전문가용 (티스토리)", "😊 환자용 (네이버)"]
        )
        target_code = "doctor" if "전문가" in target_type else "patient"

        if st.button("✍️ 글 자동 생성"):
            if b_papers.empty: st.error("논문 없음")
            elif not openai_api_key: st.error("API Key 필요")
            else:
                with st.spinner("작성 중..."):
                    article = generate_blog_article(b_date_str, b_papers, openai_api_key, b_model, target_code)
                    save_blog_post(b_date_str, target_code, article)
                    st.success("완료!")
                    st.rerun()

    with c_b2:
        st.subheader("📄 미리보기")
        tab_doc, tab_pat = st.tabs(["👨‍⚕️ 전문가용", "😊 환자용"])
        with tab_doc:
            doc_post = get_blog_post(b_date_str, "doctor")
            if doc_post:
                st.markdown(doc_post)
                st.divider()
                st.code(doc_post, language='markdown')
            else: st.info("없음")
        with tab_pat:
            pat_post = get_blog_post(b_date_str, "patient")
            if pat_post:
                st.markdown(pat_post)
                st.divider()
                st.code(pat_post, language='markdown')
            else: st.info("없음")

# --- [Tab 3: 보관함] ---
with tab_archive:
    df_all = pd.read_sql("SELECT * FROM papers", sqlite3.connect(DB_NAME))

    if df_all.empty:
        st.info("보관함이 비어있습니다. 검색 탭에서 논문을 추가해주세요.")
    else:
        st.subheader("🔍 보관함 필터링")
        
        all_categories = sorted(df_all['intervention_category'].unique().tolist())
        selected_categories = st.multiselect(
            "중재법 선택 (다중 선택 가능)", 
            all_categories, 
            default=all_categories
        )

        if 'archive_body_part' not in st.session_state:
            st.session_state.archive_body_part = "전체"

        def get_archive_btn_color(part):
            return "primary" if st.session_state.archive_body_part == part else "secondary"

        st.markdown("##### 신체 부위별 보기")
        col_b1, col_b2, col_b3, col_b4, col_b5, col_b6 = st.columns(6)
        
        parts_map = {
            "두경부": "🧠", "척추/허리": "🦴", "상지": "💪", 
            "하지": "🦵", "내장기/전신": "🫀", "전체": "🔍"
        }
        
        for i, (part, emoji) in enumerate(parts_map.items()):
            col = [col_b1, col_b2, col_b3, col_b4, col_b5, col_b6][i]
            with col:
                if st.button(f"{emoji} {part}", key=f"arc_btn_{i}", type=get_archive_btn_color(part), use_container_width=True):
                    st.session_state.archive_body_part = part
                    st.rerun()

        df_filtered = df_all.copy()
        
        if selected_categories:
            df_filtered = df_filtered[df_filtered['intervention_category'].isin(selected_categories)]
        else:
            df_filtered = pd.DataFrame() 

        if st.session_state.archive_body_part != "전체":
            df_filtered = df_filtered[df_filtered['target_body_part'] == st.session_state.archive_body_part]

        st.divider()

        st.subheader(f"📚 논문 목록 ({len(df_filtered)}건)")
        
        if not df_filtered.empty:
            df_filtered.insert(0, "delete_sel", False)
            df_filtered["url"] = "https://pubmed.ncbi.nlm.nih.gov/" + df_filtered["pmid"]
            
            edited_df = st.data_editor(
                df_filtered,
                column_config={
                    "delete_sel": st.column_config.CheckboxColumn("삭제", width="small"),
                    "url": st.column_config.LinkColumn("Link", display_text="🔗", width="small"),
                    "title_kr": st.column_config.TextColumn("제목", width="large"),
                    "summary": st.column_config.TextColumn("요약", width="large"),
                    "date_published": st.column_config.TextColumn("수집일", width="medium"),
                    "clinical_score": st.column_config.NumberColumn("점수", format="%d점", width="small"),
                    "full_text_status": st.column_config.TextColumn("분석출처", width="medium"),
                },
                column_order=["delete_sel", "url", "clinical_score", "title_kr", "intervention_category", "summary", "full_text_status", "date_published"],
                hide_index=True,
                use_container_width=True,
                key="archive_editor"
            )

            to_delete = edited_df[edited_df["delete_sel"] == True]
            if not to_delete.empty:
                st.warning(f"{len(to_delete)}건의 논문을 선택하셨습니다.")
                if st.button("🗑️ 선택한 논문 영구 삭제", type="primary"):
                    delete_papers(to_delete['pmid'].tolist())
                    st.success("삭제 완료!")
                    st.rerun()
        else:
            st.warning("조건에 맞는 논문이 없습니다.")

# --- [Tab 4: 검색] ---
with tab_search:
    col1, col2 = st.columns(2)
    with col1: s_date = st.date_input("시작", value=datetime.now()-timedelta(days=2))
    with col2: e_date = st.date_input("종료", value=datetime.now())
    
    # [수정] max_results 슬라이더 복구 및 명확한 위치 배치
    max_results = st.select_slider("검색할 최대 논문 수", options=[10, 30, 50, 100, 300], value=50)

    if 'search_res' not in st.session_state: st.session_state.search_res = None

    if st.button("1. PubMed 검색 (무료)"):
        with st.spinner("검색 중..."):
            st.session_state.search_res = search_pubmed_raw(s_date, e_date, max_results)
    
    if st.session_state.search_res:
        df_res = pd.DataFrame(st.session_state.search_res)
        df_res.insert(0, "Sel", ~df_res['is_saved'])
        
        st.subheader(f"검색 결과: {len(df_res)}건")
        edited_res = st.data_editor(df_res, column_config={"Sel": st.column_config.CheckboxColumn("선택")}, hide_index=True)
        targets = edited_res[edited_res["Sel"]]
        
        # [수정] 버튼 이름 변경
        if st.button(f"🚀 2. 선택한 {len(targets)}건 AI 분석 및 저장"):
            if not openai_api_key: st.error("Key Missing")
            else:
                conn = sqlite3.connect(DB_NAME)
                cur = conn.cursor()
                bar = st.progress(0)
                target_pmids = targets['pmid'].tolist()
                full_list = [p for p in st.session_state.search_res if p['pmid'] in target_pmids]
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
                st.success("저장 완료!")
                st.session_state.search_res = None
                time.sleep(1)
                st.rerun()
                
                
# --- [Tab 5: 자동화 설정 (NEW)] ---
with tab_settings:
    st.subheader("🤖 자동화 봇 제어판")
    st.markdown("""
    여기서 설정을 변경하면, 매일 아침 GitHub Actions 봇이 이 설정을 확인하고 작동합니다.
    (내 컴퓨터가 꺼져 있어도 돌아갑니다.)
    """)
    
    current_status = get_config("auto_bot_enabled")
    is_on = True if current_status == "True" else False
    
    st.divider()
    
    col1, col2 = st.columns([1, 2])
    with col1:
        st.markdown("#### 🚀 자동 브리핑 전송")
        auto_toggle = st.toggle("매일 아침 자동 실행 켜기", value=is_on)
        
        if auto_toggle != is_on:
            set_config("auto_bot_enabled", str(auto_toggle))
            st.success(f"설정이 저장되었습니다: {'켜짐' if auto_toggle else '꺼짐'}")
            time.sleep(1)
            st.rerun()
            
    with col2:
        st.info(f"""
        **현재 상태:** {'🟢 작동 중' if auto_toggle else '🔴 정지됨'}
        
        **작동 원리:**
        1. 매일 아침 8시(KST)에 봇이 깨어납니다.
        2. 이 DB 파일을 열어서 **'켜짐'** 상태인지 확인합니다.
        3. 켜져 있다면:
           - 어제 나온 논문을 검색합니다.
           - Top 7 논문을 심층 분석하여 **이 DB에 저장**합니다.
           - 데일리 브리핑을 작성해 **텔레그램으로 전송**합니다.
           - **업데이트된 DB를 자동으로 저장소에 백업**합니다.
        """)
        
# (나머지 메인 실행 코드)
if __name__ == "__main__":
    init_db() # 실행 시 DB 초기화/업데이트