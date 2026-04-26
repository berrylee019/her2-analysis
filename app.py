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
    """복잡한 GDC 컬럼명을 유연하게 탐색하여 HER2 데이터를 로드합니다."""
    # 형님이 올려주신 파일명 우선 탐색
    paths = [
        'data/2026-04-26T13-44_export.csv',
        '2026-04-26T13-44_export.csv',
        'data/clinical.tsv',
        'clinical.tsv'
    ]
    
    target_path = next((p for p in paths if os.path.exists(p)), None)
    if not target_path: return None
    
    try:
        # 파일 형식 자동 판별
        sep = '\t' if target_path.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_path, sep=sep)
        
        # [핵심] 컬럼명 유연 탐색 (패턴 매칭)
        cols = df_cli.columns.tolist()
        
        # IHC: 'her2'와 'ihc'가 들어간 컬럼 혹은 그냥 'her2'가 들어간 컬럼 중 가장 적절한 것 선택
        ihc_col = next((c for c in cols if 'her2' in c.lower() and 'ihc' in c.lower()), 
                       next((c for c in cols if 'her2_status' in c.lower() and 'fish' not in c.lower()), None))
        
        # FISH: 'her2'와 'fish'가 동시에 들어간 컬럼
        fish_col = next((c for c in cols if 'her2' in c.lower() and 'fish' in c.lower()), None)
        
        # ID: 'submitter_id'가 포함된 컬럼 (cases.submitter_id 등)
        id_col = next((c for c in cols if 'submitter_id' in c.lower()), None)

        def check_her2_low(row):
            # 컬럼이 없을 경우를 대비한 안전장치
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 판정 로직 (다양한 GDC 표기 대응)
            if '1+' in ihc or ihc == '1.0' or ihc == '1': return True
            if '2+' in ihc or ihc == '2.0' or ihc == '2':
                # FISH가 음성이거나 정보가 없는 경우 Low로 분류
                if any(x in fish for x in ['NEG', 'NON', 'NOT']): return True
                if fish in ['', 'NAN', '--', "'--", 'UNKNOWN', 'PND']: return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            # 매칭 실패 시 사이드바에 현재 찾은 최선의 컬럼명 노출 (디버깅용)
            st.sidebar.warning(f"탐색된 컬럼명 -> IHC: {ihc_col}, ID: {id_col}")
            return None
    except Exception as e:
        st.sidebar.error(f"파일 분석 오류: {e}")
        return None

# --- 변이 데이터 가져오기 (이전과 동일) ---
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
        data = []
        for h in hits:
            aa = h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            for occ in h.get('occurrence', []):
                case_id = occ.get('case', {}).get('submitter_id')
                data.append({"Case_ID": str(case_id).strip(), "AA_Change": aa})
        return pd.DataFrame(data)
    except: return None

# --- 에너지 분석 및 PDB 관련 함수 ---
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

# 3. 메인 화면 구성
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
her2_low_only = st.sidebar.checkbox("HER2-Low 환자군만 보기")

with st.spinner('데이터 분석 중...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 4. 데이터 통합 및 필터링
if df_mut is not None:
    df_display = df_mut.copy()
    
    if her2_low_only:
        if df_clinical is not None:
            low_patients = df_clinical[df_clinical['is_her2_low'] == True]
            low_ids = low_patients['Match_ID'].unique()
            st.sidebar.info(f"분류된 HER2-Low 환자: {len(low_ids)}명")
            
            # 실제 매칭 수행
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            if not df_display.empty:
                st.sidebar.success(f"매칭 성공: {len(df_display)}건")
            else:
                st.sidebar.warning("ID 매칭 실패 (임상-변이 ID 불일치)")
        else:
            st.sidebar.error("임상 데이터를 로드할 수 없습니다.")

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
        st.subheader("🔬 3D Structure Analysis")
        if not df_display.empty:
            mats = [m for m in df_display['AA_Change'].unique() if m != 'N/A']
            selected_mut = st.selectbox("분석할 변이 선택:", mats)
            if selected_mut:
                status, icon, dist = estimate_binding_energy(selected_mut)
                st.metric("결합 영향도", f"{icon} {status}")
                pdb_path = get_pdb_file('3WZE')
                if pdb_path: st_molstar(pdb_path, height=400)
else:
    st.error("GDC API 연결 실패")
