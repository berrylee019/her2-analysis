import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수
@st.cache_data
def load_clinical_data():
    """GDC CSV의 복잡한 컬럼 구조를 강제로 파싱하여 HER2 데이터를 추출합니다."""
    # 형님의 파일 경로
    paths = [
        'data/2026-04-26T13-44_export.csv',
        '2026-04-26T13-44_export.csv',
        'data/clinical.tsv',
        'clinical.tsv'
    ]
    
    target_path = next((p for p in paths if os.path.exists(p)), None)
    if not target_path:
        return None
    
    try:
        # CSV 로드 (헤더가 복잡할 수 있으므로 쉼표로 명확히 지정)
        df_cli = pd.read_csv(target_path, sep=None, engine='python', low_memory=False)
        
        # [핵심] 컬럼명 전수 조사 (그물망 로직)
        cols = df_cli.columns.tolist()
        
        # 1. IHC 컬럼: 이름에 'her2'와 'ihc'가 모두 들어간 모든 컬럼 찾기 (숫자/점 무관)
        ihc_candidates = [c for c in cols if 'her2' in c.lower() and 'ihc' in c.lower()]
        if not ihc_candidates:
            # ihc라는 단어가 없으면 her2_status가 포함된 거라도 검색
            ihc_candidates = [c for c in cols if 'her2_status' in c.lower() and 'fish' not in c.lower()]
        
        ihc_col = ihc_candidates[0] if ihc_candidates else None
        
        # 2. FISH 컬럼: her2와 fish가 동시에 들어간 컬럼
        fish_candidates = [c for c in cols if 'her2' in c.lower() and 'fish' in c.lower()]
        fish_col = fish_candidates[0] if fish_candidates else None
        
        # 3. ID 컬럼: submitter_id 포함 컬럼
        id_col = next((c for c in cols if 'submitter_id' in c.lower()), None)

        def check_her2_low(row):
            # 1.0, 2.0 등 숫자형과 '1+', '2+' 등 문자형 모두 대응
            ihc_val = str(row.get(ihc_col, '')).strip().upper()
            fish_val = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 판정 로직
            # IHC 1+ 계열
            if any(x in ihc_val for x in ['1+', '1.0', '1']):
                return True
            # IHC 2+ 계열 + FISH 음성
            if any(x in ihc_val for x in ['2+', '2.0', '2']):
                neg_terms = ['NEG', 'NON', 'NOT', '0']
                if any(t in fish_val for t in neg_terms): return True
                # FISH 정보가 없으면 보수적으로 Low에 포함 (GDC 특성)
                if fish_val in ['', 'NAN', '--', "'--", 'UNKNOWN', 'PND']: return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            # 여전히 못 찾을 경우 상세 리포트 (사이드바)
            st.sidebar.error(f"❌ 매칭 실패 상세: IHC={ihc_col}, ID={id_col}")
            return None
    except Exception as e:
        st.sidebar.error(f"⚠️ 파일 분석 오류: {e}")
        return None

# --- GDC 변이 데이터 호출 (기존과 동일) ---
@st.cache_data
def get_her2_mutations():
    ssm_url = "https://api.gdc.cancer.gov/ssms"
    filters = {"op": "and", "content": [
        {"op": "in", "content": {"field": "occurrence.case.project.project_id", "value": ["TCGA-BRCA"]}},
        {"op": "in", "content": {"field": "genes.symbol", "value": ["ERBB2"]}}
    ]}
    params = {"filters": json.dumps(filters), "fields": "consequence.transcript.aa_change,occurrence.case.submitter_id", "format": "JSON", "size": "1000"}
    try:
        r = requests.get(ssm_url, params=params)
        hits = r.json()['data']['hits']
        data = [{"Case_ID": str(occ.get('case', {}).get('submitter_id')).strip(), 
                 "AA_Change": h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')}
                for h in hits for occ in h.get('occurrence', [])]
        return pd.DataFrame(data)
    except: return None

# 3. 메인 UI
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
low_only = st.sidebar.checkbox("HER2-Low 환자군 분석")

with st.spinner('데이터를 매칭 중입니다...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 4. 결과 출력
if df_mut is not None:
    df_display = df_mut.copy()
    
    if low_only:
        if df_clinical is not None:
            low_ids = df_clinical[df_clinical['is_her2_low'] == True]['Match_ID'].unique()
            st.sidebar.success(f"✅ HER2-Low: {len(low_ids)}명 식별")
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            st.sidebar.info(f"🎯 매칭된 변이: {len(df_display)}건")
        else:
            st.sidebar.error("임상 데이터를 불러올 수 없습니다.")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            counts = df_display['AA_Change'].value_counts().reset_index()
            counts.columns = ['Mutation', 'Count']
            # Deprecation 경고 방지: width='stretch' 대신 use_container_width=True
            st.dataframe(counts.head(10), use_container_width=True)
        else:
            st.warning("데이터가 없습니다.")

    with col2:
        st.subheader("🔬 3D Structure Analysis")
        if not df_display.empty:
            mats = [m for m in df_display['AA_Change'].unique() if m != 'N/A']
            selected = st.selectbox("변이 선택:", mats)
            if selected:
                res_num = "".join(filter(str.isdigit, str(selected)))
                dist = abs(int(res_num) - 755) if res_num else 0
                st.metric("Dist from Pocket", f"{dist} residues")
                pdb_path = "3WZE.pdb"
                if os.path.exists(pdb_path): st_molstar(pdb_path, height=400)
                else: st.info("PDB 파일을 준비 중입니다...")
else:
    st.error("GDC API 연결 실패")
