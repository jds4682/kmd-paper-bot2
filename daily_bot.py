import os
import requests
import json
import sqlite3
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from Bio import Entrez
from openai import OpenAI
import re
import time  # [중요] 시간 지연을 위해 필수

# ===================== [환경 변수] =====================
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
EMAIL_ADDRESS = os.environ.get("EMAIL_ADDRESS")
TELEGRAM_TOKEN = os.environ.get("TELEGRAM_TOKEN")
CHAT_ID = os.environ.get("CHAT_ID")
DB_NAME = 'kmd_papers_v5_column.db' 

# 키 확인 (보안상 앞 5자리만 출력)
print(f"DEBUG: API Key Loaded? {'Yes' if OPENAI_API_KEY else 'No'}")
if OPENAI_API_KEY:
    print(f"DEBUG: Key starts with: {OPENAI_API_KEY[:5]}...")

if not TELEGRAM_TOKEN or not CHAT_ID or not OPENAI_API_KEY:
    print("❌ 설정 오류: Secrets(API KEY, Telegram Info)를 확인하세요.")
    exit(1)

Entrez.email = EMAIL_ADDRESS
client = OpenAI(api_key=OPENAI_API_KEY)

# ===================== [DB 관련 함수] =====================
def save_paper_to_db(data):
    """분석된 논문을 DB에 저장"""
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        cur.execute('''CREATE TABLE IF NOT EXISTS papers (
            pmid TEXT PRIMARY KEY, date_published TEXT, title_kr TEXT, intervention_category TEXT, 
            target_body_part TEXT, specific_point TEXT, study_design TEXT, clinical_score INTEGER,
            summary TEXT, original_title TEXT, abstract TEXT, icd_code TEXT, full_text_status TEXT
        )''')
        
        cur.execute('INSERT OR REPLACE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)', (
            data['pmid'], datetime.now().strftime('%Y-%m-%d'),
            data['title_kr'], "자동수집", 
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

# ===================== [분석 로직] =====================
def fetch_pmc_fulltext(pmid):
    try:
        # [Rate Limit 방지] Entrez도 너무 빠르면 차단당함
        time.sleep(1) 
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
    Text: {content[:8000]} 
    
    Output JSON format only: {{ "korean_title": "...", "clinical_score": 8, "summary": "...", "study_design": "...", "target_body_part": "...", "specific_point": "..." }}
    """
    try:
        # [핵심 수정] OpenAI 호출 전 2초 대기 (서버 과부하 방지)
        print("⏳ AI 분석 대기 중...")
        time.sleep(2) 
        
        res = client.chat.completions.create(
            model="gpt-4o-mini",
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0
        )
        
        raw_text = res.choices[0].message.content.strip()
        match = re.search(r'\{.*\}', raw_text, re.DOTALL)
        if match:
            data = json.loads(match.group())
            data['pmid'] = pmid
            data['original_title'] = title
            data['abstract'] = abstract
            data['source'] = status
            return data
        else:
            return None
    except Exception as e:
        # [핵심 수정] 에러 메시지를 더 자세히 출력
        print(f"🚨 OpenAI API 에러 발생: {type(e).__name__}")
        print(f"내용: {e}")
        return None

def send_telegram(msg):
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    try:
        res = requests.post(url, json={"chat_id": CHAT_ID, "text": msg, "parse_mode": "Markdown"})
        if res.status_code != 200:
            print(f"텔레그램 전송 실패: {res.text}")
    except Exception as e:
        print(f"텔레그램 연결 에러: {e}")

# ===================== [메인 실행] =====================
if __name__ == "__main__":
    print("🤖 봇 기동...")
    print("🟢 자동화 강제 실행 모드")
    
    today = datetime.now()
    yesterday = (today - timedelta(days=1)).strftime("%Y/%m/%d")
    two_days_ago = (today - timedelta(days=2)).strftime("%Y/%m/%d")
    
    term = '("TCM" OR "Acupuncture" OR "Herbal medicine") AND (hasabstract[text]) AND ("Humans"[Mesh])'
    
    print(f"🔍 검색 기간: {two_days_ago} ~ {yesterday}")
    
    try:
        handle = Entrez.esearch(db="pubmed", term=term, mindate=two_days_ago, maxdate=yesterday, datetype="pdat", retmax=7)
        pmids = Entrez.read(handle)["IdList"]
    except Exception as e:
        print(f"검색 에러: {e}")
        pmids = []
    
    if not pmids:
        print("논문 없음")
        exit(0)
        
    print(f"📄 검색된 논문 {len(pmids)}건 분석 시작...")

    analyzed_list = []
    for pmid in pmids:
        try:
            h = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
            art = Entrez.read(h)['PubmedArticle'][0]['MedlineCitation']['Article']
            title = art['ArticleTitle']
            abst = art['Abstract']['AbstractText'][0] if 'Abstract' in art else ""
            
            result = analyze_paper_bot(title, abst, pmid)
            if result:
                save_paper_to_db(result)
                analyzed_list.append(result)
            else:
                print(f"❌ {pmid} 분석 실패 (결과 없음)")
                
        except Exception as e:
            print(f"Skip {pmid}: {e}")
            
    if analyzed_list:
        analyzed_list.sort(key=lambda x: x.get('clinical_score', 0), reverse=True)
        
        briefing = f"📅 **{yesterday} 한의 임상 브리핑**\n\n"
        for i, paper in enumerate(analyzed_list[:5]):
            briefing += f"{'🥇' if i==0 else '🥈' if i==1 else '📰'} **{paper['korean_title']}**\n"
            briefing += f"(⭐{paper.get('clinical_score',0)})\n"
            briefing += f"{paper.get('summary','')}\n"
            briefing += f"🔗 https://pubmed.ncbi.nlm.nih.gov/{paper['pmid']}\n\n"
            
        send_telegram(briefing)
        print("✅ 텔레그램 전송 및 DB 저장 완료")
    else:
        print("❌ 분석된 논문이 없어 전송하지 않았습니다.")
