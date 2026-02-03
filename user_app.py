import streamlit as st
import sqlite3
import pandas as pd
from datetime import datetime
import db_handler as db  # [중요] DB 핸들러

# ===================== [앱 시작 시 DB 동기화] =====================
if 'db_synced' not in st.session_state:
    try:
        with st.spinner("최신 데이터 불러오는 중..."):
            db.pull_db()
        st.session_state.db_synced = True
    except Exception as e:
        pass

# ===================== [설정] =====================
st.set_page_config(page_title="한의학 논문 AI 큐레이터", layout="wide", page_icon="🏥")
DB_NAME = 'kmd_papers_v5_column.db'

# 관리자 로그인 정보
try:
    ADMIN_ID = st.secrets["admin"]["id"]
    ADMIN_PW = st.secrets["admin"]["pw"]
except Exception as e:
    # 로컬 테스트용 (secrets 파일이 없을 때 대비) 또는 에러 처리
    #st.error("관리자 설정(Secrets)이 없습니다.")
    
# 스타일 커스텀
st.markdown("""
    <style>
    .main-footer {text-align: center; color: grey; padding: 20px; font-size: 0.8rem;}
    .comment-box {background-color: #f0f2f6; padding: 10px; border-radius: 10px; margin-bottom: 10px;}
    .reply-box {background-color: #e8eef9; padding: 10px; border-radius: 10px; margin-left: 30px; margin-bottom: 10px; border-left: 3px solid #4e8cff;}
    </style>
""", unsafe_allow_html=True)

# ===================== [DB 구조 보정 (핵심 수정 부분)] =====================
def init_user_db():
    """테이블이 없으면 강제로 생성하는 함수"""
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    
    # 1. 댓글 테이블 생성
    cur.execute('''
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
    
    # 2. 게시판 테이블 생성
    cur.execute('''
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

# ===================== [DB 조회/저장 함수] =====================
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

# 댓글/게시판 관련
def add_comment(target_id, target_type, author, content, parent_id=None):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    # 테이블이 없을 경우를 대비해 여기서도 한 번 더 체크해도 좋음
    cur.execute("INSERT INTO comments (target_id, target_type, author, content, parent_id, created_at) VALUES (?, ?, ?, ?, ?, ?)",
                (target_id, target_type, author, content, parent_id, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
    conn.commit()
    conn.close()
    db.push_db() # GitHub 저장

def get_comments(target_id, target_type):
    conn = sqlite3.connect(DB_NAME)
    try:
        df = pd.read_sql("SELECT * FROM comments WHERE target_id=? AND target_type=? ORDER BY created_at ASC", conn, params=(target_id, target_type))
    except: 
        df = pd.DataFrame()
    conn.close()
    return df

def delete_item(table, item_id):
    conn = sqlite3.connect(DB_NAME)
    cur = conn.cursor()
    cur.execute(f"DELETE FROM {table} WHERE id=?", (item_id,))
    conn.commit()
    conn.close()
    db.push_db() # GitHub 저장

# ===================== [UI 컴포넌트] =====================

def sidebar_admin_login():
    """사이드바 관리자 로그인 창"""
    with st.sidebar:
        st.header("🔧 설정")
        if 'is_admin' not in st.session_state:
            st.session_state.is_admin = False

        if not st.session_state.is_admin:
            with st.expander("🔒 관리자 로그인"):
                uid = st.text_input("ID")
                upw = st.text_input("PW", type="password")
                if st.button("로그인"):
                    if uid == ADMIN_ID and upw == ADMIN_PW:
                        st.session_state.is_admin = True
                        st.success("관리자 모드 활성화")
                        st.rerun()
                    else:
                        st.error("정보가 올바르지 않습니다.")
        else:
            st.success("🔑 관리자님 환영합니다!")
            if st.button("로그아웃"):
                st.session_state.is_admin = False
                st.rerun()

def comment_section(target_id, target_type):
    st.subheader("💬 의견 나누기")
    
    with st.form(f"c_form_{target_id}"):
        c1, c2 = st.columns([1, 4])
        author = c1.text_input("닉네임", placeholder="예: 한의사 김")
        content = c2.text_input("내용", placeholder="자유롭게 의견을 남겨주세요.")
        if st.form_submit_button("댓글 등록"):
            if author and content:
                add_comment(target_id, target_type, author, content)
                st.rerun()
            else:
                st.warning("닉네임과 내용을 모두 입력해주세요.")

    comments = get_comments(target_id, target_type)
    if not comments.empty:
        parents = comments[comments['parent_id'].isnull()]
        for _, p in parents.iterrows():
            st.markdown(f"<div class='comment-box'><b>{p['author']}</b> <span style='color:grey;font-size:0.8em'>({p['created_at']})</span><br>{p['content']}</div>", unsafe_allow_html=True)
            
            if st.session_state.is_admin:
                if st.button("🗑️ 삭제", key=f"del_{p['id']}"):
                    delete_item('comments', p['id'])
                    st.rerun()

            with st.expander("↳ 답글 달기"):
                with st.form(f"r_form_{p['id']}"):
                    r_auth = st.text_input("닉네임", key=f"ra_{p['id']}")
                    r_cont = st.text_input("내용", key=f"rc_{p['id']}")
                    if st.form_submit_button("답글 등록"):
                        if r_auth and r_cont:
                            add_comment(target_id, target_type, r_auth, r_cont, p['id'])
                            st.rerun()
            
            children = comments[comments['parent_id'] == p['id']]
            for _, c in children.iterrows():
                st.markdown(f"<div class='reply-box'><b>↳ {c['author']}</b> <span style='color:grey;font-size:0.8em'>({c['created_at']})</span><br>{c['content']}</div>", unsafe_allow_html=True)
                if st.session_state.is_admin:
                    if st.button("🗑️ 삭제", key=f"del_{c['id']}"):
                        delete_item('comments', c['id'])
                        st.rerun()

# ===================== [메인 페이지] =====================
def main():
    # 1. 앱 시작 시 테이블 생성 (누락 방지)
    init_user_db()
    
    # 2. 관리자 로그인 처리
    sidebar_admin_login()
    
    st.title("🏥 한의학 논문 AI 큐레이터")
    st.caption("매일 업데이트되는 근거중심 한의학 정보")
    
    menu = st.tabs(["📅 데일리 브리핑", "📖 전문가 칼럼", "📚 논문 보관함", "🗣️ 자유게시판"])
    
    # --- [1] 데일리 브리핑 ---
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
            st.warning("해당 날짜의 브리핑이 아직 발행되지 않았습니다.")

    # --- [2] 전문가 칼럼 ---
    with menu[1]:
        c_date = st.date_input("발행일", datetime.now(), key="bd")
        c_str = c_date.strftime("%Y-%m-%d")
        c_type = st.radio("보기 모드", ["👨‍⚕️ 전문가용", "😊 환자용"], horizontal=True)
        target = "doctor" if "전문가" in c_type else "patient"
        
        post = get_blog_post(c_str, target)
        if post:
            st.markdown(post)
        else:
            st.info("해당 날짜의 칼럼이 없습니다.")

    # --- [3] 논문 보관함 (필터 강화) ---
    with menu[2]:
        st.subheader("📚 근거중심 한의학 논문 DB")
        df = get_papers()
        
        if not df.empty:
            with st.expander("🔎 검색 및 필터 설정", expanded=True):
                c1, c2, c3 = st.columns(3)
                search_txt = c1.text_input("제목/내용 검색")
                try:
                    min_date = pd.to_datetime(df['date_published']).min().date()
                    max_date = pd.to_datetime(df['date_published']).max().date()
                except:
                    min_date, max_date = datetime.now().date(), datetime.now().date()
                date_range = c2.date_input("연구 기간", [min_date, max_date])
                
                all_designs = sorted(df['study_design'].astype(str).unique().tolist())
                sel_designs = c3.multiselect("연구 설계 (SR, RCT 등)", all_designs)
                
                all_cats = sorted(df['intervention_category'].astype(str).unique().tolist())
                sel_cats = st.multiselect("중재법 (침, 한약 등)", all_cats)

            df_filt = df.copy()
            if search_txt:
                df_filt = df_filt[df_filt['title_kr'].str.contains(search_txt, case=False) | 
                                  df_filt['summary'].str.contains(search_txt, case=False)]
            if len(date_range) == 2:
                s_d, e_d = date_range
                df_filt['date_published'] = pd.to_datetime(df_filt['date_published']).dt.date
                df_filt = df_filt[(df_filt['date_published'] >= s_d) & (df_filt['date_published'] <= e_d)]
            if sel_designs:
                df_filt = df_filt[df_filt['study_design'].isin(sel_designs)]
            if sel_cats:
                df_filt = df_filt[df_filt['intervention_category'].isin(sel_cats)]

            st.markdown(f"**검색 결과:** 총 {len(df_filt)}건")
            st.divider()

            for _, row in df_filt.iterrows():
                with st.expander(f"[{row['study_design']}] {row['title_kr']} ({row['intervention_category']})"):
                    st.markdown(f"**임상점수:** ⭐{row['clinical_score']} | **발행일:** {row['date_published']}")
                    st.info(f"💡 **요약:** {row['summary']}")
                    st.caption(f"상세 중재: {row['specific_point']}")
                    st.link_button("PubMed 원문 보기", f"https://pubmed.ncbi.nlm.nih.gov/{row['pmid']}")
        else:
            st.info("보관된 논문이 없습니다.")

    # --- [4] 자유게시판 ---
    with menu[3]:
        st.subheader("🗣️ 자유게시판")
        st.caption("한의학 관련 자유로운 의견을 남겨주세요. (로그인 불필요)")
        
        with st.expander("📝 새 글 쓰기"):
            with st.form("board_form"):
                b_auth = st.text_input("작성자 (닉네임)")
                b_title = st.text_input("제목")
                b_content = st.text_area("내용")
                if st.form_submit_button("등록하기"):
                    if b_auth and b_title and b_content:
                        conn = sqlite3.connect(DB_NAME); cur = conn.cursor()
                        cur.execute("INSERT INTO community_board (title, content, author, created_at) VALUES (?, ?, ?, ?)", 
                                    (b_title, b_content, b_auth, datetime.now().strftime('%Y-%m-%d %H:%M:%S')))
                        conn.commit(); conn.close()
                        db.push_db() # 저장
                        st.success("등록되었습니다.")
                        st.rerun()
                    else:
                        st.warning("모든 항목을 입력해주세요.")
        
        conn = sqlite3.connect(DB_NAME)
        try:
            bdf = pd.read_sql("SELECT * FROM community_board ORDER BY created_at DESC", conn)
        except: bdf = pd.DataFrame()
        conn.close()
        
        for _, row in bdf.iterrows():
            with st.container(border=True):
                c1, c2 = st.columns([8, 1])
                with c1:
                    st.markdown(f"#### {row['title']}")
                    st.caption(f"작성자: {row['author']} | {row['created_at']}")
                    st.text(row['content'])
                with c2:
                    if st.session_state.is_admin:
                        if st.button("🗑️", key=f"bd_{row['id']}", help="게시글 삭제"):
                            delete_item('community_board', row['id'])
                            st.rerun()

    st.markdown("---")
    st.markdown('<div class="main-footer">ⓒ 2026 한의학 논문 AI 큐레이터 | 제작자: 장석우</div>', unsafe_allow_html=True)

if __name__ == "__main__":
    main()
