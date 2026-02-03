import os
import requests
import json
import sqlite3
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from Bio import Entrez
from openai import OpenAI

# ===================== [환경 변수 로드] =====================
# GitHub Actions에 등록할 비밀키들을 불러옵니다.
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
EMAIL_ADDRESS = os.environ.get("EMAIL_ADDRESS")
TELEGRAM_TOKEN = os.environ.get("TELEGRAM_TOKEN")
CHAT_ID = os.environ.get("CHAT_ID")

Entrez.email = EMAIL_ADDRESS
client = OpenAI(api_key=OPENAI_API_KEY)

# ===================== [핵심 기능 함수들] =====================
# (기존 app.py의 로직을 그대로 가져오되, Streamlit 관련 코드는 제거)

def fetch_pmc_fulltext(pmid):
    try:
        link_results = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        if not link_results or not link_results[0]['LinkSetDb']: return None, "Abstract Only"
        pmc_id = link_results[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        root = ET.fromstring(handle.read())
        full_text = "".join([t for t in root.itertext()])
        return full_text[:20000], "✅ Full Text (PMC)"
    except: return None, "Error"

def search_papers_yesterday():
    # 어제 날짜 기준 검색
    yesterday = datetime.now() - timedelta(days=1)
    date_str = yesterday.strftime("%Y/%m/%d")
    
    search_term = """
    ("TCM" OR "Traditional chinese medicine" OR "Herbal medicine" OR "Acupuncture" OR "Chuna") 
    AND (hasabstract[text]) AND ("Humans"[Mesh]) 
    AND ("Case Reports"[ptyp] OR "Clinical Trial"[ptyp] OR "Randomized Controlled Trial"[ptyp] OR "Systematic Review"[ptyp])
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, mindate=date_str, maxdate=date_str, datetype="pdat", retmax=10)
        record = Entrez.read(handle)
        return record["IdList"]
    except: return []

def analyze_and_generate_briefing(id_list):
    if not id_list: return "오늘은 새로운 임상 논문이 없습니다."
    
    analyzed_data = []
    # 상위 5개만 분석
    for pmid in id_list[:5]: 
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
            article = Entrez.read(handle)['PubmedArticle'][0]
            title = article['MedlineCitation']['Article']['ArticleTitle']
            abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText'][0]
            
            full_text, status = fetch_pmc_fulltext(pmid)
            content = full_text if full_text else abstract
            
            # PICO 분석 (GPT-4o-mini로 요약)
            response = client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{
                    "role": "user", 
                    "content": f"이 논문을 한글로 PICO 요약해줘. (약어는 풀네임으로)\nTitle: {title}\nContent: {content[:10000]}"
                }]
            )
            summary = response.choices[0].message.content
            analyzed_data.append(f"🔹 **[{title}]**\n(소스: {status})\n{summary}\n🔗 https://pubmed.ncbi.nlm.nih.gov/{pmid}\n")
        except: continue

    # 최종 브리핑 작성 (GPT-4o)
    today_str = datetime.now().strftime("%Y-%m-%d")
    prompt = f"""
    아래 논문 요약본들을 바탕으로 '한의사 {today_str} 데일리 임상 브리핑'을 작성해줘.
    텔레그램 메시지용이므로 가독성 좋게, 이모지 사용해서 작성해.
    
    [논문 데이터]
    {"".join(analyzed_data)}
    """
    
    final_res = client.chat.completions.create(
        model="gpt-4o",
        messages=[{"role": "user", "content": prompt}]
    )
    return final_res.choices[0].message.content

# ===================== [전송 함수] =====================
def send_telegram(message):
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    payload = {"chat_id": CHAT_ID, "text": message, "parse_mode": "Markdown"}
    requests.post(url, json=payload)

# ===================== [실행부] =====================
if __name__ == "__main__":
    print("🚀 자동화 봇 시작...")
    pmids = search_papers_yesterday()
    print(f"검색된 논문: {len(pmids)}건")
    
    if pmids:
        briefing_text = analyze_and_generate_briefing(pmids)
        print("브리핑 생성 완료. 전송 중...")
        send_telegram(briefing_text)
        # 여기에 post_to_blog(briefing_text) 함수만 추가하면 블로그도 자동 업로드 됨
        print("✅ 전송 완료!")
    else:
        send_telegram("오늘은 검색된 임상 논문이 없습니다. 푹 쉬세요! ☕")