import streamlit as st
import sqlite3
import pandas as pd
import json
import re
import requests
import urllib.parse
from Bio import Entrez
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from openai import OpenAI
import time
import PyPDF2
import db_handler as db

# ===================== [0. 초기 설정] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터 Pro", layout="wide", page_icon="🏥")

if 'db_synced' not in st.session_state:
    try:
        with st.spinner("데이터 동기화 중..."):
            db.pull_db()
        st.session_state.db_synced = True
    except Exception as e:
        st.warning(f"DB 동기화 실패: {e}")

# ===================== [1. 무료 파이썬 분석기] =====================
def analyze_metadata_free(text, title):
    text_lower = text.lower() if text else ""
    meta = {"n_count": "", "p_value": "", "tags": []}

    n_match = re.search(r'\bn\s*=\s*(\d+)', text_lower)
    if n_match: meta["n_count"] = f"n={n_match.group(1)}"

    if "p<0.05" in text_lower.replace(" ", "") or "p < 0.05" in text_lower:
        meta["p_value"] = "✅ P<0.05"

    if "acupotomy" in text_lower or "miniscalpel" in text_lower: meta["tags"].append("#도침")
    if "pharmacopuncture" in text_lower or "bee venom" in text_lower: meta["tags"].append("#약침")
    if "chuna" in text_lower or "tuina" in text_lower: meta["tags"].append("#추나")
    if "thread" in text_lower and "embedding" in text_lower: meta["tags"].append("#매선")
    if "herbal" in text_lower or "decoction" in text_lower: meta["tags"].append("#한약")
    
    meta["tags_str"] = ", ".join(meta["tags"])
    return meta

# ===================== [2. NEW: 추가된 심층 분석 도구들] =====================
def read_pdf_file(uploaded_file):
    try:
        reader = PyPDF2.PdfReader(uploaded_file)
        text = ""
        for page in reader.pages:
            text += page.extract_text() + "\n"
        return text[:30000]
    except: return None

# [업그레이드] 링크 생성을 위해 URL 정보도 함께 반환하도록 수정됨
def get_consensus_evidence(topic_query, required_keywords=[]):
    """
    주제(Query)로 검색하되, 결과 논문 제목에 required_keywords(핵심단어)가 
    하나라도 포함되어 있지 않으면 '가짜 결과'로 간주하고 버립니다.
    """
    try:
        # 검색어: (주제) AND (RCT/Review) AND (최신 5년)
        search_term = f"({topic_query}) AND (Systematic Review[ptyp] OR Meta-Analysis[ptyp] OR Randomized Controlled Trial[ptyp]) AND (\"2015\"[Date - Publication] : \"3000\"[Date - Publication])"
        
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=10, sort="relevance") # 10개 넉넉히 가져옴
        record = Entrez.read(handle)
        id_list = record["IdList"]
        
        if not id_list: return "관련된 추가 근거 논문이 검색되지 않았습니다.", []

        handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="xml")
        records = Entrez.read(handle)
        
        evidence_text = ""
        ref_list = []
        valid_count = 0
        
        for article in records['PubmedArticle']:
            try:
                title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_list = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                abstract = " ".join(abstract_list) if abstract_list else ""
                
                # [🛡️ 핵심: 안전장치] 제목에 핵심 키워드가 있는지 검사
                # 예: required_keywords = ['allergic', 'rhinitis']
                # 제목이 "Kidney disease..." 이면 -> 탈락!
                is_relevant = False
                if not required_keywords: # 키워드 없으면 그냥 통과
                    is_relevant = True
                else:
                    for keyword in required_keywords:
                        if keyword.lower() in title.lower():
                            is_relevant = True
                            break
                
                if not is_relevant:
                    continue # 키워드 없으면 스킵하고 다음 논문 봄

                # 검증 통과한 논문만 추가
                valid_count += 1
                evidence_text += f"\n[Ref {valid_count}] {title}\n요약: {abstract[:200]}...\n"
                
                pmid = str(article['MedlineCitation']['PMID'])
                ref_list.append({
                    "index": valid_count,
                    "title": title,
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}"
                })
                
                if valid_count >= 5: break # 5개 채우면 종료
                
            except: continue
            
        if valid_count == 0:
            return "검색 결과가 있었으나, 주제와 정확히 일치하는 논문은 없었습니다.", []
            
        return evidence_text, ref_list
        
    except Exception as e: return f"교차 검증 중 오류: {str(e)}", []

# ===================== [UI 사이드바] =====================
with st.sidebar:
    st.header("⚙️ 기본 설정")
    openai_api_key = st.text_input("OpenAI API Key", type="password")
    email_address = st.text_input("Email (PubMed용)", value="your_email@example.com")
    st.divider()
    telegram_token = st.text_input("Bot Token", type="password")
    chat_id = st.text_input("Chat ID")

Entrez.email = email_address
DB_NAME = 'kmd_papers_v5_column.db' 

# ===================== [3. DB 관리 & 마이그레이션] =====================
def migrate_db():
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    new_columns = [("n_count", "TEXT"), ("p_value", "TEXT"), ("tags", "TEXT"), ("user_note", "TEXT")]
    for col, dtype in new_columns:
        try: cursor.execute(f"ALTER TABLE papers ADD COLUMN {col} {dtype}")
        except sqlite3.OperationalError: pass
    conn.commit(); conn.close()

def get_papers_by_date(target_date_str):
    conn = sqlite3.connect(DB_NAME)
    try:
        query = "SELECT * FROM papers WHERE date_published = ?"
        df = pd.read_sql(query, conn, params=(target_date_str,))
    except: df = pd.DataFrame()
    conn.close()
    return df

def get_daily_column(date_str):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    try:
        cursor.execute("SELECT content FROM daily_columns WHERE date_id = ?", (date_str,))
        res = cursor.fetchone()
        return res[0] if res else None
    except: return None
    finally: conn.close()

def save_daily_column(date_str, content):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO daily_columns VALUES (?, ?, ?)", (date_str, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit(); conn.close(); db.push_db()

def get_blog_post(date_str, target_type):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    try:
        cursor.execute("SELECT content FROM blog_posts WHERE date_id = ? AND target_type = ?", (date_str, target_type))
        res = cursor.fetchone()
        return res[0] if res else None
    except: return None
    finally: conn.close()

def save_blog_post(date_str, target_type, content):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    cursor.execute("INSERT OR REPLACE INTO blog_posts VALUES (?, ?, ?, ?)", (date_str, target_type, content, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit(); conn.close(); db.push_db()

def delete_papers(pmid_list):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    if pmid_list:
        placeholders = ', '.join('?' for _ in pmid_list)
        cursor.execute(f"DELETE FROM papers WHERE pmid IN ({placeholders})", pmid_list)
        conn.commit()
    conn.close(); db.push_db()

def check_if_exists(pmid):
    conn = sqlite3.connect(DB_NAME); cursor = conn.cursor()
    cursor.execute("SELECT 1 FROM papers WHERE pmid=?", (pmid,))
    exists = cursor.fetchone() is not None
    conn.close()
    return exists

# ===================== [4. Full Text Fetcher] =====================
def fetch_pmc_fulltext(pmid):
    try:
        link_results = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        if not link_results or not link_results[0]['LinkSetDb']: return None, "Abstract Only"
        pmc_id = link_results[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        xml_data = handle.read()
        root = ET.fromstring(xml_data)
        full_text = "".join([text for body in root.findall(".//body") for text in body.itertext()])
        return full_text[:30000], "✅ PMC 전문(Full Text) 분석됨"
    except Exception as e: return None, f"Error: {str(e)}"

# ===================== [5. AI 분석 로직] =====================
def analyze_paper_strict(paper_data, api_key):
    client = OpenAI(api_key=api_key)
    prompt = f"""
    너는 임상 한의학 논문 분류 전문가다.
    [필수 규칙] 중재법(침/뜸/부항/한약/약침/추나/기타), 부위(두경부/척추/상지/하지/내장기)
    [JSON 형식] {{
        "korean_title": "한글 제목", "study_design": "연구 유형", "intervention_category": "카테고리",
        "target_body_part": "신체부위", "specific_point": "상세 중재 내용", "clinical_score": 8,
        "summary": "3줄 요약", "icd_code": "코드", "full_text_status": "Abstract Check"
    }}
    Title: {paper_data['title']}
    Abstract: {paper_data['abstract']}
    """
    try:
        response = client.chat.completions.create(model="gpt-4o-mini", messages=[{"role": "user", "content": prompt}], temperature=0.0)
        data = json.loads(re.search(r'\{.*\}', response.choices[0].message.content.strip(), re.DOTALL).group())
        if "DROP" in str(data.get("study_design", "")): return {"error": "DROP"}
        return data
    except Exception as e: return {"error": str(e)}

# ===================== [6. PubMed 검색] =====================
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
        id_list = Entrez.read(handle)["IdList"]
        if not id_list: return []
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
                "pmid": pmid, "title": title, "abstract": abstract,
                "predicted_category": simple_keyword_classify(title + abstract),
                "is_saved": check_if_exists(pmid)
            })
        except: continue
    return raw_papers

# ===================== [7. 데일리 브리핑 생성기] =====================
def generate_daily_briefing_pro_v3(date_str, papers_df, api_key, model_choice):
    client = OpenAI(api_key=api_key)
    top_papers = papers_df.sort_values(by='clinical_score', ascending=False).head(10)
    if top_papers.empty: return "분석할 논문이 없습니다."

    analyzed_data = []
    prog_bar = st.progress(0); status_text = st.empty()

    for idx, (_, row) in enumerate(top_papers.iterrows()):
        prog_bar.progress((idx+1)/len(top_papers))
        status_text.text(f"🔍 안전 분석 모드 동작 중... ({idx+1}): {row['title_kr']}")
        
        full_text, ft_status = fetch_pmc_fulltext(row['pmid'])
        content_source = full_text if full_text else row['abstract']
        
        pico_prompt = f"""이 논문을 PICO 구조로 분석하라. Title: {row['title_kr']} Text: {content_source[:15000]}"""
        try: pico_res_text = client.chat.completions.create(model="gpt-4o-mini", messages=[{"role": "user", "content": pico_prompt}], temperature=0.0).choices[0].message.content
        except: pico_res_text = "분석 실패"

        analyzed_data.append({
            "pmid": row['pmid'], "title": row['title_kr'], "score": row['clinical_score'],
            "study_design": row['study_design'], "source": ft_status, "detail_analysis": pico_res_text
        })
    status_text.empty(); prog_bar.empty()

    final_prompt = f"""
    [역할]
    당신은 20년 차 임상 한의사이자, 최신 의학 트렌드를 분석하는 수석 에디터입니다.
    동료 한의사들에게 매일 아침 영양가 높은 최신 논문 정보를 브리핑해주는 역할을 맡았습니다.

    [작성 목표]
    제공된 JSON 데이터를 바탕으로 가독성 좋고, 임상에 바로 적용 가능한 **'Daily Clinical Briefing'**을 작성하세요.
    
    [입력 데이터]
    {json.dumps(analyzed_data, ensure_ascii=False)}

    [작성 가이드라인]
    1. **톤앤매너:** 전문적이지만 딱딱하지 않게. "합니다/했습니다" 보다는 **"~임/함/확인됨"** 등의 개조식 혹은 **"~하는 것이 좋겠습니다"** 식의 제안형 어조 혼용.
    2. **핵심 강조:** 중요한 수치(N수, P값, 효과 크기)나 결론은 **굵은 글씨**로 강조.
    3. **비판적 시각:** 연구 디자인(RCT, SR 등)과 임상 점수(score)를 참고하여, 신뢰도가 낮은 논문은 "추가 검증이 필요함"이라고 언급.
    4. Huzhang과 같이 다른나라말에 기반한 말을 한국어로 뭔지 모르는것은 반드시 찾아서 괄호안에 주석을 달아줄것. 찾지못하였다면 오늘의 논문으로 싣는 것을 포기할 것
    5. EA, MA와 같은 줄임말을 사용한것은 EA(electroacupuncture)등과 같이 반드시 full name을 적어줄것. 찾지못하였다면 오늘의 논문으로 싣는 것을 포기할 것

    [출력 포맷 (Markdown)]
    
    # 📅 {date_str} 한의 임상 브리핑
    
    ---
    ### 🏆 Today's Pick: 집중 탐구 (Top 2)
    *가장 임상 가치가 높거나 흥미로운 논문 2개를 선정하여 상세히 다룹니다.*

    #### 1. [한글 논문 제목]
    > **연구설계:** [Study Design] | **임상점수:** ⭐[Score]/10
    
    * **🎯 핵심 요약:** (연구의 목적과 결론을 1문장으로)
    * **💊 중재/방법:** (무엇을 어떻게 했는지? 구체적인 혈자리, 처방명, 용량 등)
    * **📊 주요 결과:** (대조군 대비 구체적인 수치 변화 포함)
    * **💡 임상 제언:** (이 결과를 실제 한의원 진료실에서 어떻게 해석하고 적용해야 할지 에디터의 관점으로 서술. 예: "만성 요통 환자 티칭 시 참고할 만합니다.")
    * 🔗 [원문 보기](https://pubmed.ncbi.nlm.nih.gov/[PMID])

    (2번 논문도 동일 형식)

    ---
    ### 📰 Research Shorts: 놓치면 아쉬운 단신 (5선)
    *제목과 핵심 결과만 빠르게 훑어봅니다.*

    1. **[한글 제목]** ([Study Design])
       - 결과: 핵심 결과 3줄 요약
       - 🔗 [Link](https://pubmed.ncbi.nlm.nih.gov/[PMID])
    
    2. (나머지 논문들...)

    ---
   (Data source: PubMed)_
    """
    try: return client.chat.completions.create(model=model_choice, messages=[{"role": "user", "content": final_prompt}]).choices[0].message.content
    except Exception as e: return f"생성 실패: {e}"

# ===================== [9. 메인 UI] =====================
migrate_db() # 앱 실행 시 DB 구조 자동 업데이트

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
                st.success("완료!")
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
                if not telegram_token or not chat_id: st.error("토큰 필요")
                else:
                    try:
                        url = f"https://api.telegram.org/bot{telegram_token}/sendMessage"
                        res = requests.post(url, json={"chat_id": chat_id, "text": final_msg, "parse_mode": "Markdown"})
                        if res.status_code == 200: st.success("전송 완료")
                        else: st.error(f"실패: {res.text}")
                    except Exception as e: st.error(f"에러: {e}")
            st.divider()
            st.markdown(final_msg)
        else: st.warning("브리핑 없음")
# --- [Tab 2: 블로그 (참고문헌 자동 추가 기능)] ---
# --- [Tab 2: 블로그 (상위 개념 확장 검색 기능 탑재)] ---
with tab_blog:
    c_b1, c_b2 = st.columns([1, 3])
    with c_b1:
        st.subheader("✒️ 심층 블로그 생성")
        b_date = st.date_input("날짜", value=datetime.now(), key="blog_date")
        b_date_str = b_date.strftime("%Y-%m-%d")
        b_papers = get_papers_by_date(b_date_str)
        st.info(f"후보: {len(b_papers)}건")
        
        if not b_papers.empty:
            sel_title = st.selectbox("논문 선택", b_papers['title_kr'].tolist())
            target_paper = b_papers[b_papers['title_kr'] == sel_title].iloc[0]
            uploaded_pdf = st.file_uploader("📄 PDF 업로드 (선택)", type="pdf")
            b_model = st.selectbox("모델:", ["gpt-4o", "gpt-4o-mini"], index=0)
            target_type = st.radio("타겟:", ["👨‍⚕️ 전문가용", "😊 환자용"])

            if st.button("🚀 심층 분석 & 글쓰기"):
                if not openai_api_key: st.error("Key 없음")
                else:
                    # 1. 본문 텍스트 확보
                    with st.spinner("1. 자료 분석 중... (PDF/PMC)"):
                        status_msg = "초록(Abstract) 기반"
                        content_source = target_paper['abstract']
                        if uploaded_pdf:
                            pdf_txt = read_pdf_file(uploaded_pdf)
                            if pdf_txt: content_source = pdf_txt; status_msg = "📂 PDF 전문 분석"
                        elif not uploaded_pdf:
                            pmc_txt, pmc_msg = fetch_pmc_fulltext(target_paper['pmid'])
                            if pmc_txt: content_source = pmc_txt; status_msg = pmc_msg

                    # 2. [핵심 변경] 검색어 2개 생성 (Specific / Broad)
                    with st.spinner("2. 확장형 교차 검증(Consensus) 실행 중..."):
                        client = OpenAI(api_key=openai_api_key)
                        
                        # AI에게 "좁은 검색어"와 "넓은 검색어"를 동시에 요청
                        q_prompt = f"""
                        Analyze this text and generate TWO English search queries for PubMed validation.
                        
                        1. "Specific": Exact intervention + Exact disease (e.g., "Sopoongsan AND Atopic Dermatitis")
                        2. "Broad": Intervention class + Symptom/Disease category (e.g., "Herbal Medicine AND Pruritus", "TCM AND Eczema")
                        
                        * Only output valid JSON.
                        
                        Text: {content_source[:1500]}
                        
                        Output JSON format:
                        {{
                            "specific": "Query string 1",
                            "broad": "Query string 2"
                        }}
                        """
                        try:
                            q_resp = client.chat.completions.create(model="gpt-4o-mini", messages=[{"role":"user","content":q_prompt}]).choices[0].message.content
                            q_json = json.loads(re.search(r'\{.*\}', q_resp, re.DOTALL).group())
                            
                            q_specific = q_json['specific']
                            q_broad = q_json['broad']
                            
                            # 1차 검색 (좁은 범위)
                            evidence, ref_list = get_consensus_evidence(q_specific)
                            
                            # [로직] 결과가 3개 미만이면 -> 넓은 범위 검색 추가 실행!
                            if len(ref_list) < 3:
                                ev_broad, ref_broad = get_consensus_evidence(q_broad)
                                evidence += f"\n\n[추가 근거 (상위/유사 계열)]: {q_broad}\n" + ev_broad
                                
                                # 중복 제거하며 리스트 합치기
                                existing_urls = [r['url'] for r in ref_list]
                                for item in ref_broad:
                                    if item['url'] not in existing_urls:
                                        # 인덱스 번호 조정
                                        item['index'] = len(ref_list) + 1
                                        ref_list.append(item)
                                        
                        except Exception as e:
                            evidence = f"검증 데이터 생성 실패: {e}"
                            ref_list = []

                    # 3. 글 작성 (비교/보완 관점 추가)
                    with st.spinner("3. 심층 칼럼 작성 중..."):
                        final_prompt = f"""
                        당신은 임상 한의학 전문 작가입니다. 
                        제공된 [논문]과 [검증자료]를 종합하여 블로그 글을 작성하세요.

                        [상황] {status_msg}
                        [메인 논문] {content_source[:25000]}
                        
                        [검증 자료 (Consensus)] 
                        {evidence}

                        [작성 핵심 지침 - 비교와 보완]
                        [Role & Persona] 당신은 20년 경력의 베테랑 임상 한의사이자, 난해한 의학 논문을 일반인의 언어로 쉽고 재치 있게 풀어내는 '건강 전문 칼럼니스트'입니다. 신뢰감 있는 말투를 유지하되, 네이버 블로그 특유의 친근함과 가독성을 갖춘 글을 작성하세요.

[Context & Task] 제기된 주제에서의 학술적 내용을 바탕으로, 한의사가 실질적인 도움을 받을 수 있는 블로그 포스팅을 작성해 주세요.

[Specific Guidelines]

타겟 독자: 메인 중재에 대한 질환이 가장 잘 발생하는 연령대 및 성별을 대상으로 메인 중재의 대상질환 환자 및 그들의 자녀 혹은 가족 (가독성이 높고 문장이 명확해야 함).

구조(Framework): '문제 제기 -> 새로운 대안 제시(연구 결과) -> 핵심 약재 설명(한의학 중재) -> 전문가의 비평 -> 결론' 순으로 구성하세요.

데이터 보완: > * 제공된 레퍼런스 외에 핵심 중재의 **Network Pharmacology(네트워크 약리학)**적 기전(예: 염증성 사이토카인 억제 등) 또는 국가제공 데이터, 제공되지 않았지만 확실하게 확인된 임상데이터 혹은 논문을 조사해서 한개 ~두개의 문장으로 쉽게 덧붙여주세요. 매번 네트워크 약리학만 근거로 적을 수는 없으니 다른 근거를 넣는것도 항상 고려해주세요.

만약에 논문에서 지적된 '낮은 방법론적 질'을 언급할 때, 독자가 불안해하지 않도록 앞으로 더 정교한 연구가 기대되는 유망한 분야라는 점을 강조하세요.

톤앤매너: 신뢰(80%) + 다정함(20%). 이모지를 적절히 사용하고, 전문 용어는 반드시 괄호나 쉬운 비유로 풀어서 설명하세요.

SEO 최적화: 메인 치료법, 메인 치료법의 대상 질환 치료, 해당 중재의 효능 등의 키워드가 자연스럽게 녹아들게 하세요.
                        5. 마지막에는 [Reference] 리스트를 달지 마세요. (시스템이 알아서 답니다.)
                        """
                        article = client.chat.completions.create(model=b_model, messages=[{"role":"user","content":final_prompt}]).choices[0].message.content
                        
                        # 참고문헌 부착
                        if ref_list:
                            article += "\n\n---\n### 📚 참고 문헌 (References)\n"
                            article += f"**[Main]** {target_paper['title_kr']} (PMID: {target_paper['pmid']})\n"
                            for ref in ref_list:
                                article += f"{ref['index']}. [{ref['title']}]({ref['url']})\n"

                        save_blog_post(b_date_str, "doctor" if "전문가" in target_type else "patient", article)
                        st.success(f"작성 완료! (참고문헌 {len(ref_list)}건 확보)")
                        st.rerun()

    with c_b2:
        st.subheader("📄 미리보기")
        t1, t2 = st.tabs(["👨‍⚕️ 전문가용", "😊 환자용"])
        with t1:
            post = get_blog_post(b_date_str, "doctor")
            if post: st.markdown(post)
        with t2:
            post = get_blog_post(b_date_str, "patient")
            if post: st.markdown(post)

# --- [Tab 3: 보관함] ---
with tab_archive:
    df_all = pd.read_sql("SELECT * FROM papers", sqlite3.connect(DB_NAME))
    if not df_all.empty:
        st.subheader("🔍 필터링")
        cats = sorted(df_all['intervention_category'].unique().tolist())
        sel_cats = st.multiselect("중재법", cats, default=cats)
        df_filt = df_all[df_all['intervention_category'].isin(sel_cats)] if sel_cats else df_all.copy()
        
        if not df_filt.empty:
            df_filt.insert(0, "del", False)
            df_filt["url"] = "https://pubmed.ncbi.nlm.nih.gov/" + df_filt["pmid"]
            
            edited = st.data_editor(
                df_filt,
                column_config={
                    "del": st.column_config.CheckboxColumn("삭제", width="small"),
                    "url": st.column_config.LinkColumn("Link"),
                    "tags": st.column_config.TextColumn("태그"),
                    "n_count": st.column_config.TextColumn("N수"),
                    "user_note": st.column_config.TextColumn("메모")
                },
                column_order=["del", "url", "clinical_score", "tags", "n_count", "title_kr", "user_note"],
                hide_index=True, use_container_width=True
            )
            
            if st.button("🗑️ 삭제 확인"):
                to_del = edited[edited["del"]]['pmid'].tolist()
                if to_del: delete_papers(to_del); st.rerun()

            st.divider()
            st.caption("👇 논문 상세 검증")
            for _, row in df_filt.iterrows():
                with st.expander(f"{row['title_kr']}"):
                    st.info(row['summary'])
                    c1, c2, c3 = st.columns(3)
                    c1.link_button("📄 원문", row['url'], use_container_width=True)
                    q = urllib.parse.quote(row['original_title'])
                    c2.link_button("⚖️ Consensus", f"https://consensus.app/results/?q={q}", use_container_width=True)
                    c3.link_button("🤖 SciSpace", f"https://typeset.io/search?q={q}", use_container_width=True)

# --- [Tab 4: 검색] ---
with tab_search:
    c1, c2 = st.columns(2)
    s_d = c1.date_input("시작", datetime.now()-timedelta(days=2))
    e_d = c2.date_input("끝", datetime.now())
    limit = st.slider("개수", 10, 100, 50)
    
    if 'search_res' not in st.session_state: st.session_state.search_res = None
    if st.button("1. 검색"):
        with st.spinner(".."): st.session_state.search_res = search_pubmed_raw(s_d, e_d, limit)
        
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
                
                # 테이블 생성
                cur.execute('''CREATE TABLE IF NOT EXISTS papers (
                    pmid TEXT PRIMARY KEY, date_published TEXT, title_kr TEXT, 
                    intervention_category TEXT, target_body_part TEXT, specific_point TEXT, 
                    study_design TEXT, clinical_score INTEGER, summary TEXT, 
                    original_title TEXT, abstract TEXT, icd_code TEXT, full_text_status TEXT,
                    n_count TEXT, p_value TEXT, tags TEXT, user_note TEXT
                )''')
                
                full_list = [p for p in st.session_state.search_res if p['pmid'] in targets['pmid'].tolist()]
                for i, p in enumerate(full_list):
                    bar.progress((i+1)/len(full_list))
                    res = analyze_paper_strict(p, openai_api_key)
                    meta = analyze_metadata_free(p['abstract'], p['title'])
                    
                    if "error" not in res:
                        cur.execute("INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                            p['pmid'], datetime.now().strftime('%Y-%m-%d'),
                            res.get('korean_title'), res.get('intervention_category'),
                            res.get('target_body_part'), res.get('specific_point'),
                            res.get('study_design'), res.get('clinical_score'),
                            res.get('summary'), p['title'], p['abstract'], 
                            res.get('icd_code'), "Abstract Saved",
                            meta['n_count'], meta['p_value'], meta['tags_str'], ""
                        ))
                        conn.commit()
                conn.close(); db.push_db()
                st.success("완료!"); st.session_state.search_res = None; st.rerun()

if __name__ == "__main__":
    if not st.session_state.get('db_synced'):
        db.pull_db()
        st.session_state.db_synced = True
    migrate_db()










