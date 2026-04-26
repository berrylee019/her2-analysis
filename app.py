import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수 정의
@st.cache_data
def load_clinical_data():
    """모든 경로와 다양한 컬럼명을 자동으로 탐색하여 임상 데이터를 로드합니다."""
    paths = [
        'data/2026-04-26T13-44_export.csv',
        '2026-04-26T13-44_export.csv',
        'data/clinical.tsv',
        'clinical.tsv'
    ]
    
    target_path = next((p for p in paths if os.path.exists(p)), None)
    if not target_path: return None
    
    try:
        # 파일 형식 자동 판별 및 로드
        sep = '\t' if target_path.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_path, sep=sep, low_memory=False)
        
        # [핵심] 그물망 컬럼 탐색
        cols = df_cli.columns.tolist()
        
        # IHC 컬럼: 'her2'와 'ihc'가 모두 들어간 컬럼을 찾되, 없으면 'her2_status' 포함 컬럼 탐색
        ihc_col = next((c for c in cols if 'her2' in c.lower() and 'ihc' in c.lower()), 
                       next((c for c in cols if 'her2_status' in c.lower() and 'fish' not in c.lower()), None))
        
        # FISH 컬럼: 'her2'와 'fish'가 모두 들어간 컬럼
        fish_col = next((c for c in cols if 'her2' in c.lower() and 'fish' in c.lower()), None)
        
        # ID 컬럼: 'submitter_id'가 포함된 컬럼 (보통 cases.submitter_id)
        id_col = next((c for c in cols if 'submitter_id' in c.lower()), 
                      next((c for c in cols if 'case_id' in c.lower()), None))

        def check_her2_low(row):
            # 데이터를 문자열로 변환 후 공백 제거 및 대문자화
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 판정 (숫자형 1.0, 2.0 및 텍스트 1+, 2+ 모두 대응)
            # 1) IHC 1인 경우
            if any(x in ihc for x in ['1+', '1.0', '1']): return True
            
            # 2) IHC 2인 경우 (FISH 음성 확인)
            if any(x in ihc for x in ['2+', '2.0', '2']):
                # FISH가 음성이거나 정보가 없는 경우(GDC 특성상 미기입이 많음)를 Low로 판단
                neg_keywords = ['NEG', 'NON', 'NOT', '0']
                if any(k in fish for k in neg_keywords): return True
                if fish in ['', 'NAN', '--', "'--", 'UNKNOWN', 'PND', 'NOT PERFORMED']: return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            # 실패 시 사이드바에 디버깅 정보 출력
            st.sidebar.warning(f"⚠️ 컬럼 탐색 실패\n- 찾은 IHC컬럼: {ihc_col}\n- 찾은 ID컬럼: {id_col}")
            return None
    except Exception as e:
        st.sidebar.error(f"파일 처리 오류: {e}")
        return None

# --- 이하 API 및 분석 로직 (기존과 동일하되 최적화) ---
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

def estimate_binding_energy(res_num_str):
    try:
        res_int = int("".join(filter(str.isdigit, str(res_num_str))))
        dist = abs(res_int - 755)
        if dist < 15: return "Critical", "🔴", dist
        elif dist < 30: return "Moderate", "🟡", dist
        else: return "Low", "🟢", dist
    except: return "Unknown", "⚪", 0

def get_pdb_file(pdb_id):
    file_path = f"{pdb_id}.pdb"
    if not os.path.exists(file_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        r = requests.get(url)
        with open(file_path, "w") as f: f.write(r.text)
    return file_path

# 3. 메인 UI
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
her2_low_only = st.sidebar.checkbox("HER2-Low 환자군만 보기", value=False)

with st.spinner('데이터를 불러오는 중...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 4. 필터링 및 결과 출력
if df_mut is not None:
    df_display = df_mut.copy()
    
    if her2_low_only:
        if df_clinical is not None:
            low_ids = df_clinical[df_clinical['is_her2_low'] == True]['Match_ID'].unique()
            st.sidebar.info(f"✅ HER2-Low 환자군: {len(low_ids)}명 확보")
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            if not df_display.empty:
                st.sidebar.success(f"🎯 매칭된 변이 데이터: {len(df_display)}건")
            else:
                st.sidebar.warning("매칭된 데이터가 없습니다. (ID 형식 확인 필요)")
        else:
            st.sidebar.error("임상 데이터를 로드하지 못했습니다.")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            top_mats = df_display['AA_Change'].value_counts().reset_index()
            top_mats.columns = ['Mutation', 'Count']
            st.dataframe(top_mats.head(10), use_container_width=True)
        else:
            st.warning("데이터가 없습니다.")

    with col2:
        st.subheader("🔬 3D Structure & Energy Analysis")
        if not df_display.empty:
            mats = [m for m in df_display['AA_Change'].unique() if m != 'N/A']
            selected_mut = st.selectbox("분석할 변이 선택:", mats)
            if selected_mut:
                status, icon, dist = estimate_binding_energy(selected_mut)
                st.metric("Binding Impact", f"{icon} {status} (Dist: {dist})")
                pdb_path = get_pdb_file('3WZE')
                if pdb_path: st_molstar(pdb_path, height=400)
else:
    st.error("GDC API 연결 실패")
