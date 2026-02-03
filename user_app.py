import streamlit as st
import sqlite3
import pandas as pd
import hashlib
from datetime import datetime

# ===================== [설정] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터", layout="wide", page_icon="🏥")
DB_NAME = 'kmd_papers_v5_column.db'

# 스타일 커스텀 (반응형 & 제작자 표시)
st.markdown("""
    <style>
    .main-footer {text-align: center; color: grey; padding: 20px; font-size: 0.8rem;}
    .comment-box {background-color: #f0f2f6; padding: 10px; border-radius: 10px; margin-bottom: 10px;}
    .reply-box {background-color: #e8eef9; padding: 10px; border-radius: 10px; margin-left: 30px; margin-bottom: 10px; border-left: 3px solid #4e8cff;}
    </style>
""", unsafe_allow_html=True)

# ===================== [DB 관리] =====================
def init_user_db():
    conn = sqlite3.connect(DB_NAME)
    cursor = conn.cursor()
    
    # 1. 유저 테이블
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS users (
            username TEXT PRIMARY KEY,
            password TEXT,
            nickname TEXT,
            role TEXT DEFAULT 'user'
        )
    ''')
    
    # 2. 댓글 테이블 (Target: 날짜ID 또는 게시글ID)
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS comments (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            target_id TEXT,
            target_type TEXT, 
            author TEXT,
            content TEXT,
            parent_id INTEGER,
            created_at TEXT
        )
    ''')
    
    # 3. 자유게시판 테이블
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS community_board (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            title TEXT,
            content TEXT,
            author TEXT,
            created_at TEXT
        )
    ''')
    
    conn.commit()
    conn.close()

# ===================== [인증 함수] =====================
def hash_pw(password):
    return hashlib.sha256(password.encode()).hexdigest()

def login_user(username, password):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute("SELECT nickname, role FROM users WHERE username=? AND password=?", (username, hash_pw(password)))
    user = cur.fetchone()
    conn.close()
    return user

def register_user(username, password, nickname):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        # 첫 번째 가입자는 자동으로 관리자(admin) 부여
        cur.execute("SELECT count(*) FROM users")
        cnt = cur.fetchone()[0]
        role = 'admin' if cnt == 0 else 'user'
        
        cur.execute("INSERT INTO users VALUES (?, ?, ?, ?)", (username, hash_pw(password), nickname, role))
        conn.commit()
        return True
    except:
        return False
    finally:
        conn.close()

# ===================== [데이터 조회 함수] =====================
def get_daily_briefing(date_str):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    # daily_columns 테이블이 있는지 확인 (Admin 앱에서 생성됨)
    try:
        cur.execute("SELECT content FROM daily_columns WHERE date_id = ?", (date_str,))
        res = cur.fetchone()
        return res[0] if res else None
    except: return None
    finally: conn.close()

def get_blog_post(date_str, target_type):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        cur.execute("SELECT content FROM blog_posts WHERE date_id = ? AND target_type = ?", (date_str, target_type))
        res = cur.fetchone()
        return res[0] if res else None
    except: return None
    finally: conn.close()

def get_papers():
    conn = sqlite3.connect(DB_NAME)
    try:
        df = pd.read_sql("SELECT * FROM papers", conn)
    except: df = pd.DataFrame()
    conn.close()
    return df

# ===================== [댓글/게시판 시스템] =====================
def add_comment(target_id, target_type, author, content, parent_id=None):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute("INSERT INTO comments (target_id, target_type, author, content, parent_id, created_at) VALUES (?, ?, ?, ?, ?, ?)",
                (target_id, target_type, author, content, parent_id, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()

def get_comments(target_id, target_type):
    conn = sqlite3.connect(DB_NAME)
    df = pd.read_sql("SELECT * FROM comments WHERE target_id=? AND target_type=? ORDER BY created_at ASC", conn, params=(target_id, target_type))
    conn.close()
    return df

def delete_item(table, item_id):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute(f"DELETE FROM {table} WHERE id=?", (item_id,))
    conn.commit()
    conn.close()

# ===================== [UI 컴포넌트] =====================
def sidebar_login():
    with st.sidebar:
        st.title("👤 로그인")
        if 'user' not in st.session_state:
            tab1, tab2 = st.tabs(["로그인", "회원가입"])
            with tab1:
                uid = st.text_input("아이디", key="l_id")
                upw = st.text_input("비밀번호", type="password", key="l_pw")
                if st.button("로그인"):
                    user_info = login_user(uid, upw)
                    if user_info:
                        st.session_state['user'] = {'id': uid, 'name': user_info[0], 'role': user_info[1]}
                        st.success(f"환영합니다, {user_info[0]}님!")
                        st.rerun()
                    else:
                        st.error("아이디 또는 비밀번호 오류")
            with tab2:
                new_id = st.text_input("새 아이디", key="r_id")
                new_pw = st.text_input("새 비밀번호", type="password", key="r_pw")
                new_nick = st.text_input("닉네임 (실명 X)", key="r_nick")
                if st.button("가입하기"):
                    if register_user(new_id, new_pw, new_nick):
                        st.success("가입 완료! 로그인해주세요.")
                    else:
                        st.error("이미 존재하는 아이디입니다.")
        else:
            u = st.session_state['user']
            st.info(f"👋 {u['name']}님 ({'관리자' if u['role']=='admin' else '회원'})")
            if st.button("로그아웃"):
                del st.session_state['user']
                st.rerun()

def comment_section(target_id, target_type):
    st.subheader("💬 의견 나누기")
    
    # 댓글 입력
    if 'user' in st.session_state:
        with st.form(f"c_form_{target_id}"):
            txt = st.text_area("내용을 입력하세요", height=70)
            if st.form_submit_button("댓글 달기"):
                add_comment(target_id, target_type, st.session_state['user']['name'], txt)
                st.rerun()
    else:
        st.info("로그인 후 댓글을 남길 수 있습니다.")

    # 댓글 표시
    comments = get_comments(target_id, target_type)
    if not comments.empty:
        # 부모 댓글과 자식 댓글 분리
        parents = comments[comments['parent_id'].isnull()]
        
        for _, p in parents.iterrows():
            # 부모 댓글 렌더링
            st.markdown(f"""
            <div class="comment-box">
                <b>{p['author']}</b> <span style="font-size:0.8em; color:grey;">{p['created_at']}</span><br>
                {p['content']}
            </div>
            """, unsafe_allow_html=True)
            
            # 관리자 삭제 버튼
            if 'user' in st.session_state and st.session_state['user']['role'] == 'admin':
                if st.button("🗑️ 삭제", key=f"del_{p['id']}"):
                    delete_item('comments', p['id'])
                    st.rerun()

            # 대댓글 달기 (토글)
            if 'user' in st.session_state:
                with st.expander("↳ 답글 달기"):
                    with st.form(f"r_form_{p['id']}"):
                        reply_txt = st.text_input("답글 내용")
                        if st.form_submit_button("등록"):
                            add_comment(target_id, target_type, st.session_state['user']['name'], reply_txt, p['id'])
                            st.rerun()

            # 자식 댓글 렌더링
            children = comments[comments['parent_id'] == p['id']]
            for _, c in children.iterrows():
                st.markdown(f"""
                <div class="reply-box">
                    <b>↳ {c['author']}</b> <span style="font-size:0.8em; color:grey;">{c['created_at']}</span><br>
                    {c['content']}
                </div>
                """, unsafe_allow_html=True)
                if 'user' in st.session_state and st.session_state['user']['role'] == 'admin':
                    if st.button("🗑️ 삭제", key=f"del_{c['id']}"):
                        delete_item('comments', c['id'])
                        st.rerun()

# ===================== [메인 페이지 로직] =====================
def main():
    init_user_db()
    sidebar_login()
    
    st.title("🏥 한의학 논문 AI 큐레이터")
    
    # 상단 메뉴
    menu = st.tabs(["📅 데일리 브리핑", "📖 전문가 칼럼", "📚 논문 보관함", "🗣️ 자유게시판"])
    
    # 1. 데일리 브리핑
    with menu[0]:
        sel_date = st.date_input("날짜 선택", datetime.now())
        date_str = sel_date.strftime("%Y-%m-%d")
        
        briefing = get_daily_briefing(date_str)
        if briefing:
            st.markdown("---")
            st.markdown(briefing)
            st.divider()
            comment_section(date_str, "briefing")
        else:
            st.warning("해당 날짜의 브리핑이 없습니다.")

    # 2. 전문가 칼럼 (블로그)
    with menu[1]:
        col_date = st.date_input("발행일 선택", datetime.now(), key="blog_d")
        col_str = col_date.strftime("%Y-%m-%d")
        
        col_type = st.radio("보기 모드", ["👨‍⚕️ 전문가용", "😊 환자용"], horizontal=True)
        type_code = "doctor" if "전문가" in col_type else "patient"
        
        post = get_blog_post(col_str, type_code)
        if post:
            st.markdown(post)
        else:
            st.info("발행된 칼럼이 없습니다.")

    # 3. 논문 보관함
    with menu[2]:
        st.subheader("📚 근거중심 한의학 논문 DB")
        df = get_papers()
        if not df.empty:
            # 검색 필터
            search_txt = st.text_input("논문 제목 또는 내용 검색")
            if search_txt:
                df = df[df['title_kr'].str.contains(search_txt) | df['summary'].str.contains(search_txt)]
            
            # 카테고리 필터
            cat = st.selectbox("중재법", ["전체"] + list(df['intervention_category'].unique()))
            if cat != "전체":
                df = df[df['intervention_category'] == cat]

            for _, row in df.iterrows():
                with st.expander(f"[{row['intervention_category']}] {row['title_kr']}"):
                    st.markdown(f"**임상점수:** ⭐{row['clinical_score']} | **연구유형:** {row['study_design']}")
                    st.info(row['summary'])
                    st.caption(f"발행일: {row['date_published']}")
                    st.link_button("원문 보기 (PubMed)", f"https://pubmed.ncbi.nlm.nih.gov/{row['pmid']}")
        else:
            st.info("보관된 논문이 없습니다.")

    # 4. 자유게시판
    with menu[3]:
        st.subheader("🗣️ 자유게시판")
        
        # 글쓰기
        if 'user' in st.session_state:
            with st.expander("📝 새 글 쓰기"):
                with st.form("board_form"):
                    b_title = st.text_input("제목")
                    b_content = st.text_area("내용")
                    if st.form_submit_button("등록"):
                        conn = sqlite3.connect(DB_NAME)
                        cur = conn.cursor()
                        cur.execute("INSERT INTO community_board (title, content, author, created_at) VALUES (?, ?, ?, ?)",
                                    (b_title, b_content, st.session_state['user']['name'], datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
                        conn.commit()
                        conn.close()
                        st.rerun()
        
        # 게시글 목록
        conn = sqlite3.connect(DB_NAME)
        board_df = pd.read_sql("SELECT * FROM community_board ORDER BY created_at DESC", conn)
        conn.close()
        
        for _, row in board_df.iterrows():
            with st.container(border=True):
                st.markdown(f"**{row['title']}**")
                st.caption(f"작성자: {row['author']} | {row['created_at']}")
                st.text(row['content'])
                
                # 삭제 (본인 또는 관리자)
                if 'user' in st.session_state:
                    if st.session_state['user']['role'] == 'admin' or st.session_state['user']['name'] == row['author']:
                        if st.button("삭제", key=f"bd_{row['id']}"):
                            delete_item('community_board', row['id'])
                            st.rerun()

    # Footer
    st.markdown("---")
    st.markdown('<div class="main-footer">ⓒ 2026 한의학 논문 AI 큐레이터 | 제작자: 장석우</div>', unsafe_allow_html=True)

if __name__ == "__main__":
    main()
