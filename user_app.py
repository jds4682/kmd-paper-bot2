import streamlit as st
import sqlite3
import pandas as pd
import hashlib
from datetime import datetime
import db_handler as db  # [중요] DB 핸들러 임포트

# ===================== [앱 시작 시 DB 동기화] =====================
if 'db_synced' not in st.session_state:
    with st.spinner("최신 데이터 불러오는 중..."):
        db.pull_db()
    st.session_state.db_synced = True

# ===================== [설정] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터", layout="wide", page_icon="🏥")
DB_NAME = 'kmd_papers_v5_column.db'

# 스타일 커스텀
st.markdown("""
    <style>
    .main-footer {text-align: center; color: grey; padding: 20px; font-size: 0.8rem;}
    .comment-box {background-color: #f0f2f6; padding: 10px; border-radius: 10px; margin-bottom: 10px;}
    .reply-box {background-color: #e8eef9; padding: 10px; border-radius: 10px; margin-left: 30px; margin-bottom: 10px; border-left: 3px solid #4e8cff;}
    </style>
""", unsafe_allow_html=True)

# ===================== [DB 관리] =====================
def init_user_db():
    pass # db_handler가 처리함

# 인증 관련
def hash_pw(password):
    return hashlib.sha256(password.encode()).hexdigest()

def login_user(username, password):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        cur.execute("SELECT nickname, role FROM users WHERE username=? AND password=?", (username, hash_pw(password)))
        user = cur.fetchone()
        return user
    except: return None
    finally: conn.close()

def register_user(username, password, nickname):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    try:
        cur.execute("SELECT count(*) FROM users")
        cnt = cur.fetchone()[0]
        role = 'admin' if cnt == 0 else 'user'
        
        cur.execute("INSERT INTO users VALUES (?, ?, ?, ?)", (username, hash_pw(password), nickname, role))
        conn.commit()
    except:
        return False
    finally:
        conn.close()
    
    # [중요] 회원가입 후 GitHub 업로드
    db.push_db()
    return True

# 데이터 조회
def get_daily_briefing(date_str):
    conn = sqlite3.connect(DB_NAME)
    try:
        cur = conn.cursor()
        cur.execute("SELECT content FROM daily_columns WHERE date_id = ?", (date_str,))
        res = cur.fetchone()
        return res[0] if res else None
    except: return None
    finally: conn.close()

def get_blog_post(date_str, target_type):
    conn = sqlite3.connect(DB_NAME)
    try:
        cur = conn.cursor()
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

# 댓글/게시판
def add_comment(target_id, target_type, author, content, parent_id=None):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute("INSERT INTO comments (target_id, target_type, author, content, parent_id, created_at) VALUES (?, ?, ?, ?, ?, ?)",
                (target_id, target_type, author, content, parent_id, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()
    # [중요] 댓글 저장 후 GitHub 업로드
    db.push_db()

def get_comments(target_id, target_type):
    conn = sqlite3.connect(DB_NAME)
    try:
        df = pd.read_sql("SELECT * FROM comments WHERE target_id=? AND target_type=? ORDER BY created_at ASC", conn, params=(target_id, target_type))
    except: df = pd.DataFrame()
    conn.close()
    return df

def delete_item(table, item_id):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute(f"DELETE FROM {table} WHERE id=?", (item_id,))
    conn.commit()
    conn.close()
    # [중요] 삭제 후 GitHub 업로드
    db.push_db()

# ===================== [UI] =====================
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
    if 'user' in st.session_state:
        with st.form(f"c_form_{target_id}"):
            txt = st.text_area("내용을 입력하세요", height=70)
            if st.form_submit_button("댓글 달기"):
                add_comment(target_id, target_type, st.session_state['user']['name'], txt)
                st.rerun()
    else:
        st.info("로그인 후 댓글을 남길 수 있습니다.")

    comments = get_comments(target_id, target_type)
    if not comments.empty:
        parents = comments[comments['parent_id'].isnull()]
        for _, p in parents.iterrows():
            st.markdown(f"<div class='comment-box'><b>{p['author']}</b> <span style='color:grey;font-size:0.8em'>{p['created_at']}</span><br>{p['content']}</div>", unsafe_allow_html=True)
            
            if 'user' in st.session_state and st.session_state['user']['role'] == 'admin':
                if st.button("🗑️ 삭제", key=f"del_{p['id']}"): delete_item('comments', p['id']); st.rerun()

            if 'user' in st.session_state:
                with st.expander("↳ 답글"):
                    with st.form(f"r_form_{p['id']}"):
                        rtxt = st.text_input("답글")
                        if st.form_submit_button("등록"):
                            add_comment(target_id, target_type, st.session_state['user']['name'], rtxt, p['id'])
                            st.rerun()

            children = comments[comments['parent_id'] == p['id']]
            for _, c in children.iterrows():
                st.markdown(f"<div class='reply-box'><b>↳ {c['author']}</b> <span style='color:grey;font-size:0.8em'>{c['created_at']}</span><br>{c['content']}</div>", unsafe_allow_html=True)
                if 'user' in st.session_state and st.session_state['user']['role'] == 'admin':
                    if st.button("🗑️ 삭제", key=f"del_{c['id']}"): delete_item('comments', c['id']); st.rerun()

# ===================== [메인] =====================
def main():
    init_user_db()
    sidebar_login()
    
    st.title("🏥 한의학 논문 AI 큐레이터")
    menu = st.tabs(["📅 데일리 브리핑", "📖 전문가 칼럼", "📚 논문 보관함", "🗣️ 자유게시판"])
    
    with menu[0]:
        sel_date = st.date_input("날짜", datetime.now())
        date_str = sel_date.strftime("%Y-%m-%d")
        briefing = get_daily_briefing(date_str)
        if briefing:
            st.markdown("---"); st.markdown(briefing); st.divider()
            comment_section(date_str, "briefing")
        else: st.warning("브리핑이 없습니다.")

    with menu[1]:
        c_date = st.date_input("날짜", datetime.now(), key="bd")
        c_str = c_date.strftime("%Y-%m-%d")
        c_type = st.radio("보기", ["👨‍⚕️ 전문가용", "😊 환자용"], horizontal=True)
        post = get_blog_post(c_str, "doctor" if "전문가" in c_type else "patient")
        if post: st.markdown(post)
        else: st.info("칼럼이 없습니다.")

    with menu[2]:
        st.subheader("📚 논문 DB")
        df = get_papers()
        if not df.empty:
            txt = st.text_input("검색")
            if txt: df = df[df['title_kr'].str.contains(txt) | df['summary'].str.contains(txt)]
            cat = st.selectbox("중재", ["전체"] + list(df['intervention_category'].unique()))
            if cat != "전체": df = df[df['intervention_category'] == cat]
            
            for _, row in df.iterrows():
                with st.expander(f"[{row['intervention_category']}] {row['title_kr']}"):
                    st.markdown(f"**점수:** ⭐{row['clinical_score']} | **유형:** {row['study_design']}")
                    st.info(row['summary'])
                    st.link_button("원문 보기", f"https://pubmed.ncbi.nlm.nih.gov/{row['pmid']}")
        else: st.info("보관된 논문이 없습니다.")

    with menu[3]:
        st.subheader("🗣️ 자유게시판")
        if 'user' in st.session_state:
            with st.expander("📝 글쓰기"):
                with st.form("b_form"):
                    bt = st.text_input("제목")
                    bc = st.text_area("내용")
                    if st.form_submit_button("등록"):
                        conn = sqlite3.connect(DB_NAME); cur = conn.cursor()
                        cur.execute("INSERT INTO community_board (title, content, author, created_at) VALUES (?, ?, ?, ?)", (bt, bc, st.session_state['user']['name'], datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
                        conn.commit(); conn.close()
                        # [중요] 게시글 작성 후 GitHub 업로드
                        db.push_db()
                        st.rerun()
        
        conn = sqlite3.connect(DB_NAME)
        try:
            bdf = pd.read_sql("SELECT * FROM community_board ORDER BY created_at DESC", conn)
        except: bdf = pd.DataFrame()
        conn.close()
        
        for _, row in bdf.iterrows():
            with st.container(border=True):
                st.markdown(f"**{row['title']}**")
                st.caption(f"{row['author']} | {row['created_at']}")
                st.text(row['content'])
                if 'user' in st.session_state:
                    if st.session_state['user']['role'] == 'admin' or st.session_state['user']['name'] == row['author']:
                        if st.button("삭제", key=f"bd_{row['id']}"):
                            delete_item('community_board', row['id'])
                            st.rerun()

    st.markdown("---")
    st.markdown('<div class="main-footer">ⓒ 2026 한의학 논문 AI 큐레이터 | 제작자: 장석우</div>', unsafe_allow_html=True)

if __name__ == "__main__":
    main()
