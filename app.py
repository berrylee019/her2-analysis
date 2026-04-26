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
    """형님이 업로드하신 특정 CSV 파일을 직접 로드하고 HER2-Low 환자군을 분류합니다."""
    # 형님이 알려주신 파일명과 경로를 직접 타겟팅합니다.
    # 파일이 data 폴더 안에 있는지, 혹은 실행 경로에 있는지 확인 필요
    possible_paths = [
        'data/2026-04-26T13-44_export.csv',
        '2026-04-26T13-44_export.csv'
    ]
    
    file_path = None
    for path in possible_paths:
        if os.path.exists(path):
            file_path = path
            break

    if not file_path:
        # 파일이 없을 경우 현재 디렉토리의 모든 csv 중 'export'가 포함된 파일을 찾음
        files = [f for f in os.listdir('.') if 'export.csv' in f]
        if files:
            file_path = files[0]
        else:
            return None
    
    try:
        # CSV 파일 로드
        df_cli = pd.read_csv(file_path)
        
        # [핵심] HER2 관련 컬럼 탐색 (GDC의 복잡한 컬럼명 대응)
        ihc_col = next((c for c in df_cli.columns if 'her2_status_by_ihc' in c.lower()), None)
        fish_col = next((c for c in df_cli.columns if 'her2_fish_status' in c.lower()), None)
        id_col = next((c for c in df_cli.columns if 'cases.submitter_id' in c.lower()), 
                      next((c for c in df_cli.columns if 'submitter_id' in c.lower()), None))

        def check_her2_low(row):
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 정의: IHC 1+ 또는 (IHC 2+ 이면서 FISH Negative 계열)
            if '1+' in ihc: return True
            if '2+' in ihc:
                if any(x in fish for x in ['NEGATIVE', 'NON-AMPLIFIED', 'NOT AMPLIFIED']):
                    return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            return None
    except Exception as e:
        st.error(f"데이터 로딩 중 오류: {e}")
        return None

def estimate_binding_energy(res_num_str, drug_pocket_center=755):
    try:
        res_int = int("".join(filter(str.isdigit, str(res_num_str))))
        distance = abs(res_int - drug_pocket_center)
        if distance < 15: return "Critical", "🔴", distance
        elif distance < 30: return "Moderate", "🟡", distance
        else: return "Low", "🟢", distance
    except: return "Unknown", "⚪", 0

@st.cache_data
def get_her2_mutations():
    ssm_url = "https://api.gdc.cancer.gov/ssms"
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "occurrence.case.project.project_id", "value": ["TCGA-BRCA"]}},
            {"op": "in", "content": {"field": "genes.symbol", "value": ["ERBB2"]}}
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "consequence.transcript.aa_change,occurrence.case.submitter_id",
        "format": "JSON",
        "size": "1000"
    }
    try:
        r = requests.get(ssm_url, params=params)
        hits = r.json()['data']['hits']
        data = []
        for h in hits:
            aa = h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            occurrences = h.get('occurrence', [])
            for occ in occurrences:
                case_id = occ.get('case', {}).get('submitter_id')
                data.append({"Case_ID": str(case_id).strip(), "AA_Change": aa})
        return pd.DataFrame(data)
    except: return None

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

with st.spinner('데이터를 분석 중입니다...'):
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
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            if not df_display.empty:
                st.sidebar.success(f"매칭된 변이 데이터: {len(df_display)}건")
            else:
                st.sidebar.warning("매칭된 변이 데이터가 없습니다.")
        else:
            st.sidebar.error("임상 데이터 로드 실패 (파일명을 확인하세요)")

    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            top_mats = df_display['AA_Change'].value_counts().reset_index()
            top_mats.columns = ['Mutation', 'Count']
            # 경고 해결: use_container_width=True -> width='stretch'
            st.dataframe(top_mats.head(10), width='stretch')
        else:
            st.warning("표시할 데이터가 없습니다.")
            top_mats = pd.DataFrame(columns=['Mutation', 'Count'])

    with col2:
        st.subheader("🔬 3D Structure & Energy Analysis")
        if not top_mats.empty:
            mutation_options = [m for m in top_mats['Mutation'].unique() if m != 'N/A']
            selected_mut = st.selectbox("분석할 변이를 선택하세요:", mutation_options)
            
            if selected_mut:
                res_num = "".join(filter(str.isdigit, str(selected_mut)))
                status, icon, dist = estimate_binding_energy(res_num)
                
                c1, c2 = st.columns(2)
                c1.metric("결합 영향도", f"{icon} {status}")
                c2.metric("포켓과의 거리", f"{dist} residues")
                
                pdb_path = get_pdb_file('3WZE')
                if pdb_path:
                    st_molstar(pdb_path, key='her2_viewer', height=400)
                    st.caption(f"📍 분석 지점: {selected_mut} (활성 부위 755번 기준)")
else:
    st.error("GDC API 연결 실패")
