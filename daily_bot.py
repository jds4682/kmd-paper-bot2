import os
import requests
import json
import sqlite3
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from Bio import Entrez
from openai import OpenAI

# ===================== [환경 변수] =====================
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
EMAIL_ADDRESS = os.environ.get("EMAIL_ADDRESS")
TELEGRAM_TOKEN = os.environ.get("TELEGRAM_TOKEN")
CHAT_ID = os.environ.get("CHAT_ID")
DB_NAME = 'kmd_papers_v5_column.db' 

if not TELEGRAM_TOKEN or not CHAT_ID:
    print("❌ 설정 오류: Secrets를 확인하세요.")
    exit(1)

Entrez.email = EMAIL_ADDRESS
client = OpenAI(api_key=OPENAI_API_KEY)

# ===================== [DB 관련 함수] =====================
def get_config_status():
    """Streamlit에서 설정한 자동화 ON/OFF 값을 읽어옴"""
    try:
        conn = sqlite3.connect(DB_NAME)
        cur = conn.cursor()
        cur.execute("SELECT value FROM system_config WHERE key='auto_bot_enabled'")
        res = cur.fetchone()
        conn.close()
        return res[0] == "True" if res else False
    except:
        return False # 테이블이 없거나 에러나면 안 돌림

def save_paper_to_db(data):
    """분석된 논문을 DB에 저장 (토큰 절약 핵심)"""
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        cur.execute('INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', (
            data['pmid'], datetime.now().strftime('%Y-%m-%d'),
            data['title_kr'], "자동수집", # 카테고리는 자동
            data.get('target_body_part', '기타'), data.get('specific_point', ''),
            data.get('study_design', ''), data.get('clinical_score', 0),
            data.get('summary', ''), data['original_title'], data['abstract'], 
            data.get('icd_code', ''), data.get('source', '')
        ))
        conn.commit()
        print(f"💾 DB 저장 완료: {data['title_kr']}")
    except Exception as e:
        print(f"❌ DB 저장 실패: {e}")
    finally:
        conn.close()

# ===================== [분석 로직 (app.py와 동일)] =====================
def fetch_pmc_fulltext(pmid):
    try:
        link = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        if not link or not link[0]['LinkSetDb']: return None, "Abstract Only"
        pmc_id = link[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        root = ET.fromstring(handle.read())
        return "".join([t for t in root.itertext()])[:20000], "✅ Full Text (PMC)"
    except: return None, "Error"

def analyze_paper_bot(title, abstract, pmid):
    full_text, status = fetch_pmc_fulltext(pmid)
    content = full_text if full_text else abstract
    
    prompt = f"""
    이 논문을 분석해서 JSON으로 반환해.
    [규칙]
    1. clinical_score: 1~10 (임상가치)
    2. summary: PICO 요약 (3줄)
    3. korean_title: 한글 제목 번역
    4. study_design: RCT, SR, etc.
    
    Title: {title}
    Text: {content[:10000]}
    
    Output JSON format only: {{ "korean_title": "...", "clinical_score": 8, "summary": "...", "study_design": "...", "target_body_part": "...", "specific_point": "..." }}
    """
    try:
        res = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0
        )
        data = json.loads(re.search(r'\{.*\}', res.choices[0].message.content.strip(), re.DOTALL).group())
        data['pmid'] = pmid
        data['original_title'] = title
        data['abstract'] = abstract
        data['source'] = status
        return data
    except Exception as e:
        print(f"분석 에러: {e}")
        return None

def send_telegram(msg):
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    requests.post(url, json={"chat_id": CHAT_ID, "text": msg, "parse_mode": "Markdown"})

# ===================== [메인 실행] =====================
if __name__ == "__main__":
    print("🤖 봇 기동...")
    
    # 1. 설정 확인
    if not get_config_status():
        print("🔕 자동화 설정이 꺼져있어 종료합니다.")
        exit(0)
        
    print("🟢 자동화 설정 ON - 작업 시작")
    
    # 2. 논문 검색 (어제 날짜)
    yesterday = (datetime.now() - timedelta(days=1)).strftime("%Y/%m/%d")
    term = '("TCM" OR "Acupuncture" OR "Herbal medicine") AND (hasabstract[text]) AND ("Humans"[Mesh])'
    
    try:
        handle = Entrez.esearch(db="pubmed", term=term, mindate=yesterday, maxdate=yesterday, datetype="pdat", retmax=7)
        pmids = Entrez.read(handle)["IdList"]
    except: pmids = []
    
    if not pmids:
        send_telegram(f"📅 {yesterday}\n새로운 임상 논문이 없습니다.")
        exit(0)
        
    # 3. 분석 및 DB 저장
    analyzed_list = []
    for pmid in pmids:
        try:
            h = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
            art = Entrez.read(h)['PubmedArticle'][0]['MedlineCitation']['Article']
            title = art['ArticleTitle']
            abst = art['Abstract']['AbstractText'][0] if 'Abstract' in art else ""
            
            # AI 분석
            result = analyze_paper_bot(title, abst, pmid)
            if result:
                save_paper_to_db(result) # DB에 영구 저장!
                analyzed_list.append(result)
        except Exception as e:
            print(f"Skip {pmid}: {e}")
            
    # 4. 브리핑 생성 및 전송
    if analyzed_list:
        # 점수순 정렬
        analyzed_list.sort(key=lambda x: x['clinical_score'], reverse=True)
        
        briefing = f"📅 **{yesterday} 한의 임상 브리핑**\n\n"
        for i, paper in enumerate(analyzed_list[:5]): # Top 5만
            briefing += f"{'🥇' if i==0 else '🥈' if i==1 else '📰'} **{paper['korean_title']}**\n"
            briefing += f"(⭐{paper['clinical_score']} / {paper['study_design']})\n"
            briefing += f"{paper['summary']}\n"
            briefing += f"🔗 https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}\n\n"
            
        send_telegram(briefing)
        print("✅ 전송 및 저장 완료")
