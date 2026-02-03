import streamlit as st
import sqlite3
from github import Github
import os

# 설정 가져오기
GITHUB_TOKEN = st.secrets["general"]["GH_TOKEN"]
REPO_NAME = st.secrets["general"]["REPO_NAME"]
DB_FILE = 'kmd_papers_v5_column.db'

def get_repo():
    g = Github(GITHUB_TOKEN)
    return g.get_repo(REPO_NAME)

def pull_db():
    """앱 시작할 때 GitHub에서 최신 DB 파일을 다운로드"""
    if os.path.exists(DB_FILE):
        return # 이미 있으면 패스 (속도 향상)
        
    try:
        repo = get_repo()
        contents = repo.get_contents(DB_FILE)
        with open(DB_FILE, 'wb') as f:
            f.write(contents.decoded_content)
        print("✅ GitHub에서 DB 다운로드 완료")
    except:
        print("⚠️ DB 파일이 아직 없습니다. 새로 생성합니다.")
        init_local_db()

def push_db():
    """데이터가 변경되면 GitHub에 DB 파일을 업로드 (저장)"""
    try:
        repo = get_repo()
        with open(DB_FILE, 'rb') as f:
            content = f.read()
        
        # 파일이 이미 있으면 업데이트, 없으면 생성
        try:
            contents = repo.get_contents(DB_FILE)
            repo.update_file(contents.path, "🤖 Auto-save DB", content, contents.sha)
            print("✅ DB 업데이트 완료")
        except:
            repo.create_file(DB_FILE, "🎉 Init DB", content)
            print("✅ DB 새로 생성 완료")
            
    except Exception as e:
        st.error(f"GitHub 저장 실패: {e}")

def init_local_db():
    """로컬 DB 테이블 생성 (빈 깡통 만들기)"""
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    # 논문 테이블
    c.execute('''CREATE TABLE IF NOT EXISTS papers (pmid TEXT PRIMARY KEY, date_published TEXT, title_kr TEXT, intervention_category TEXT, target_body_part TEXT, specific_point TEXT, study_design TEXT, clinical_score INTEGER, summary TEXT, original_title TEXT, abstract TEXT, icd_code TEXT, full_text_status TEXT)''')
    # 브리핑 테이블
    c.execute('''CREATE TABLE IF NOT EXISTS daily_columns (date_id TEXT PRIMARY KEY, content TEXT, created_at TEXT)''')
    # 블로그 테이블
    c.execute('''CREATE TABLE IF NOT EXISTS blog_posts (date_id TEXT, target_type TEXT, content TEXT, created_at TEXT, PRIMARY KEY (date_id, target_type))''')
    # 유저 테이블
    c.execute('''CREATE TABLE IF NOT EXISTS users (username TEXT PRIMARY KEY, password TEXT, nickname TEXT, role TEXT)''')
    # 게시판/댓글
    c.execute('''CREATE TABLE IF NOT EXISTS community_board (id INTEGER PRIMARY KEY AUTOINCREMENT, title TEXT, content TEXT, author TEXT, created_at TEXT)''')
    c.execute('''CREATE TABLE IF NOT EXISTS comments (id INTEGER PRIMARY KEY AUTOINCREMENT, target_id TEXT, target_type TEXT, author TEXT, content TEXT, parent_id INTEGER, created_at TEXT)''')
    conn.commit()
    conn.close()

# === 실행 래퍼 함수 (app.py에서 쓸 것들) ===
def run_query(query, params=(), commit=False):
    """쿼리 실행 후 저장이 필요하면 GitHub Push까지 수행"""
    # 1. 최신 DB 확인
    if not os.path.exists(DB_FILE): pull_db()
    
    conn = sqlite3.connect(DB_FILE)
    cur = conn.cursor()
    res = None
    try:
        cur.execute(query, params)
        if commit:
            conn.commit()
            push_db() # <--- 여기가 핵심! 저장하면 GitHub로 쏘기
        else:
            res = cur.fetchall()
    except Exception as e:
        print(f"DB Error: {e}")
    finally:
        conn.close()
    return res
