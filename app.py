import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정 (최상단 고정)
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수 정의
@st.cache_data
def load_clinical_data():
    """GDC CSV의 복잡한 계층형 컬럼명을 뚫고 HER2 데이터를 추출합니다."""
    # 형님이 알려주신 파일명 우선 탐색
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
        # CSV/TSV 자동 구분 로드
        sep = '\t' if target_path.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_path, sep=sep, low_memory=False)
        
        # [핵심 로직] 컬럼명 전수 조사 (패턴 매칭)
        all_cols = df_cli.columns.tolist()
        
        # 1. IHC 컬럼 찾기 (her2와 ihc가 동시에 포함된 모든 컬럼 검색)
        ihc_candidates = [c for c in all_cols if 'her2' in c.lower() and 'ihc' in c.lower()]
        # 만약 위 조건으로 못 찾으면 her2_status 포함된 거라도 검색
        if not ihc_candidates:
            ihc_candidates = [c for c in all_cols if 'her2_status' in c.lower() and 'fish' not in c.lower()]
        
        ihc_col = ihc_candidates[0] if ihc_candidates else None
        
        # 2. FISH 컬럼 찾기
        fish_candidates = [c for c in all_cols if 'her2' in c.lower() and 'fish' in c.lower()]
        fish_col = fish_candidates[0] if fish_candidates else None
        
        # 3. ID 컬럼 찾기 (cases.submitter_id 확인됨)
        id_col = next((c for c in all_cols if 'submitter_id' in c.lower()), None)

        def check_her2_low(row):
            # 값을 문자열로 변환하여 1.0, 2.0 등 숫자형 대응
            ihc_val = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish_val = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 기준: IHC 1+ 또는 (IHC 2+ 이면서 FISH 음성)
            # '1'이 포함되거나 1.0인 경우
            if any(x in ihc_val for x in ['1+', '1.0', '1']):
                return True
            # '2'가 포함되거나 2.0인 경우 + FISH 음성/미기입
            if any(x in ihc_val for x in ['2+', '2.0', '2']):
                neg_terms = ['NEG', 'NON', 'NOT', '0']
                if any(t in fish_val for t in neg_terms): return True
                if fish_val in ['', 'NAN', '--', "'--", 'PND', 'UNKNOWN']: return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            # 실패 시 사이드바에 현재 컬럼 상태 출력 (디버깅용)
            st.sidebar.warning(f"⚠️ 컬럼 탐색 결과\n- IHC: {ihc_col}\n- ID: {id_col}")
            return None
    except Exception as e:
        st.sidebar.error(f"파일 분석 오류: {e}")
        return None

# --- 변이 데이터 가져오기 (GDC API) ---
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

# --- 분석 및 시각화 관련 함수 ---
def estimate_impact(mut_str):
    try:
        res_num = int("".join(filter(str.isdigit, str(mut_str))))
        dist = abs(res_num - 755)
        if dist < 15: return "Critical", "🔴", dist
        elif dist < 30: return "Moderate", "🟡", dist
        else: return "Low", "🟢", dist
    except: return "Unknown", "⚪", 0

def get_pdb(pdb_id):
    path = f"{pdb_id}.pdb"
    if not os.path.exists(path):
        r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        with open(path, "w") as f: f.write(r.text)
    return path

# 3. 메인 UI 구성
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
low_only = st.sidebar.checkbox("HER2-Low 환자군 분석", value=False)

with st.spinner('데이터 분석 중...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 4. 데이터 통합 및 출력
if df_mut is not None:
    df_display = df_mut.copy()
    
    if low_only:
        if df_clinical is not None:
            low_patients = df_clinical[df_clinical['is_her2_low'] == True]
            low_ids = low_patients['Match_ID'].unique()
            st.sidebar.info(f"✅ HER2-Low 환자: {len(low_ids)}명 확보")
            
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            if not df_display.empty:
                st.sidebar.success(f"🎯 매칭된 변이 데이터: {len(df_display)}건")
            else:
                st.sidebar.warning("매칭된 변이 데이터가 없습니다.")
        else:
            st.sidebar.error("임상 데이터를 로드할 수 없습니다.")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            counts = df_display['AA_Change'].value_counts().reset_index()
            counts.columns = ['Mutation', 'Count']
            st.dataframe(counts.head(10), use_container_width=True)
        else:
            st.warning("데이터가 없습니다.")

    with col2:
        st.subheader("🔬 3D Structure Analysis")
        if not df_display.empty:
            mats = [m for m in df_display['AA_Change'].unique() if m != 'N/A']
            selected = st.selectbox("분석 변이:", mats)
            if selected:
                status, icon, dist = estimate_impact(selected)
                st.metric("Impact on Binding Pocket", f"{icon} {status} (Dist: {dist})")
                pdb_path = get_pdb('3WZE')
                if pdb_path: st_molstar(pdb_path, height=400)
else:
    st.error("GDC API 연결 실패")
