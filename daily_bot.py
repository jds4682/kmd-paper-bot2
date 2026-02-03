import os
import requests
import json
import xml.etree.ElementTree as ET
from datetime import datetime, timedelta
from Bio import Entrez
from openai import OpenAI

# ===================== [환경 변수 로드 & 디버깅] =====================
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY")
EMAIL_ADDRESS = os.environ.get("EMAIL_ADDRESS")
TELEGRAM_TOKEN = os.environ.get("TELEGRAM_TOKEN")
CHAT_ID = os.environ.get("CHAT_ID")

print(f"DEBUG: Token check: {'OK' if TELEGRAM_TOKEN else 'MISSING'}")
print(f"DEBUG: Chat ID check: {'OK' if CHAT_ID else 'MISSING'} (ID: {CHAT_ID})")

if not TELEGRAM_TOKEN or not CHAT_ID:
    print("❌ 에러: 텔레그램 토큰이나 채팅 ID가 설정되지 않았습니다. Settings > Secrets를 확인하세요.")
    exit(1)

Entrez.email = EMAIL_ADDRESS
client = OpenAI(api_key=OPENAI_API_KEY)

# ===================== [핵심 기능] =====================
def fetch_pmc_fulltext(pmid):
    try:
        link_results = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        if not link_results or not link_results[0]['LinkSetDb']: return None, "Abstract Only"
        pmc_id = link_results[0]['LinkSetDb'][0]['Link'][0]['Id']
        handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="xml")
        root = ET.fromstring(handle.read())
        full_text = "".join([t for t in root.itertext()])
        return full_text[:20000], "✅ Full Text (PMC)"
    except Exception as e: return None, f"Error: {e}"

def search_papers_recent():
    # [수정] 테스트를 위해 검색 기간을 '최근 3일'로 늘림
    today = datetime.now()
    start_date = today - timedelta(days=3)
    
    str_start = start_date.strftime("%Y/%m/%d")
    str_end = today.strftime("%Y/%m/%d")
    
    print(f"🔍 검색 기간: {str_start} ~ {str_end}")
    
    search_term = """
    ("TCM" OR "Traditional chinese medicine" OR "Herbal medicine" OR "Acupuncture" OR "Chuna") 
    AND (hasabstract[text]) AND ("Humans"[Mesh]) 
    AND ("Case Reports"[ptyp] OR "Clinical Trial"[ptyp] OR "Randomized Controlled Trial"[ptyp] OR "Systematic Review"[ptyp])
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=search_term, mindate=str_start, maxdate=str_end, datetype="pdat", retmax=5)
        record = Entrez.read(handle)
        return record["IdList"]
    except Exception as e:
        print(f"❌ 검색 중 에러 발생: {e}")
        return []

def analyze_and_generate_briefing(id_list):
    if not id_list: return "검색된 논문이 없습니다."
    
    analyzed_data = []
    print(f"📝 총 {len(id_list)}개 논문 분석 시작...")
    
    for pmid in id_list:
        try:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="xml")
            article = Entrez.read(handle)['PubmedArticle'][0]
            title = article['MedlineCitation']['Article']['ArticleTitle']
            abstract_list = article['MedlineCitation']['Article']['Abstract']['AbstractText']
            abstract = " ".join(abstract_list) if abstract_list else ""
            
            full_text, status = fetch_pmc_fulltext(pmid)
            content = full_text if full_text else abstract
            
            print(f" - 분석중: {title[:30]}... ({status})")
            
            response = client.chat.completions.create(
                model="gpt-4o-mini",
                messages=[{
                    "role": "user", 
                    "content": f"이 한의학 논문을 한글로 3줄 요약해줘. (PICO형식)\nTitle: {title}\nContent: {content[:5000]}"
                }]
            )
            summary = response.choices[0].message.content
            analyzed_data.append(f"🔹 **[{title}]**\n{summary}\n🔗 https://pubmed.ncbi.nlm.nih.gov/{pmid}\n")
        except Exception as e:
            print(f"   ⚠️ 개별 논문 분석 에러: {e}")
            continue

    if not analyzed_data: return "분석 중 오류가 발생하여 요약본이 없습니다."

    print("🤖 최종 브리핑 생성 중...")
    today_str = datetime.now().strftime("%Y-%m-%d")
    prompt = f"""
    아래 내용을 바탕으로 '한의사 {today_str} 데일리 임상 브리핑'을 작성해.
    
    [논문 데이터]
    {"".join(analyzed_data)}
    """
    
    final_res = client.chat.completions.create(
        model="gpt-4o",
        messages=[{"role": "user", "content": prompt}]
    )
    return final_res.choices[0].message.content

# ===================== [수정된 전송 함수 (디버깅)] =====================
def send_telegram(message):
    print("🚀 텔레그램 전송 시도 중...")
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    payload = {"chat_id": CHAT_ID, "text": message, "parse_mode": "Markdown"}
    
    try:
        response = requests.post(url, json=payload)
        # 응답 코드 확인
        if response.status_code == 200:
            print("✅ 텔레그램 전송 성공! 핸드폰을 확인하세요.")
        else:
            print(f"❌ 텔레그램 전송 실패! 상태 코드: {response.status_code}")
            print(f"❌ 에러 메시지: {response.text}")
    except Exception as e:
        print(f"❌ 연결 에러: {e}")

# ===================== [실행부] =====================
if __name__ == "__main__":
    print("🚀 자동화 봇 V2 시작...")
    pmids = search_papers_recent()
    print(f"검색된 논문 ID: {pmids}")
    
    if pmids:
        briefing_text = analyze_and_generate_briefing(pmids)
        send_telegram(briefing_text)
    else:
        print("검색된 논문이 없어 알림 메시지를 보냅니다.")
        send_telegram("오늘은 검색된 임상 논문이 없습니다. (No papers found in last 3 days)")
