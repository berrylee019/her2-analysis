import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정 (최상단)
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수 정의
@st.cache_data
def load_clinical_data():
    """폴더 내의 모든 CSV/TSV 파일을 뒤져서 임상 데이터를 자동으로 찾아냅니다."""
    # 검색할 경로들
    search_dirs = ['./', './data/']
    target_file = None
    
    # 1단계: 형님이 업로드한 특정 파일명 우선 탐색
    specific_filename = '2026-04-26T13-44_export.csv'
    for d in search_dirs:
        path = os.path.join(d, specific_filename)
        if os.path.exists(path):
            target_file = path
            break
            
    # 2단계: 못 찾으면 폴더 내의 아무 CSV나 TSV 파일을 후보로 등록
    if not target_file:
        for d in search_dirs:
            if os.path.exists(d):
                files = [f for f in os.listdir(d) if f.endswith(('.csv', '.tsv'))]
                if files:
                    target_file = os.path.join(d, files[0])
                    break

    if not target_file:
        return None
    
    try:
        # 파일 확장자에 따른 로드 (쉼표 또는 탭 구분)
        sep = '\t' if target_file.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_file, sep=sep)
        
        # [중요] HER2 관련 컬럼명 자동 매칭 (복잡한 GDC 명칭 대응)
        ihc_col = next((c for c in df_cli.columns if 'her2_status_by_ihc' in c.lower()), None)
        fish_col = next((c for c in df_cli.columns if 'her2_fish_status' in c.lower()), None)
        id_col = next((c for c in df_cli.columns if 'cases.submitter_id' in c.lower() or 'case_submitter_id' in c.lower()), 
                      next((c for c in df_cli.columns if 'submitter_id' in c.lower()), None))

        def check_her2_low(row):
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            # HER2-Low: IHC 1+ 혹은 (IHC 2+ 이면서 FISH 음성)
            if '1+' in ihc: return True
            if '2+' in ihc:
                if any(x in fish for x in ['NEGATIVE', 'NON-AMPLIFIED', 'NOT AMPLIFIED']) or fish in ['', 'NAN', 'PND']:
                    return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        return None
    except:
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
            st.sidebar.error("임상 데이터 파일을 찾을 수 없습니다. (CSV/TSV 파일을 data 폴더에 넣어주세요)")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            top_mats = df_display['AA_Change'].value_counts().reset_index()
            top_mats.columns = ['Mutation', 'Count']
            # 경고 수정 완료
            st.dataframe(top_mats.head(10), width=None, use_container_width=True)
        else:
            st.warning("데이터가 없습니다.")
            top_mats = pd.DataFrame(columns=['Mutation', 'Count'])

    with col2:
        st.subheader("🔬 3D Structure & Energy Analysis")
        if not top_mats.empty:
            mutation_options = [m for m in top_mats['Mutation'].unique() if m != 'N/A']
            selected_mut = st.selectbox("변이 선택:", mutation_options)
            if selected_mut:
                res_num = "".join(filter(str.isdigit, str(selected_mut)))
                status, icon, dist = estimate_binding_energy(res_num)
                c1, c2 = st.columns(2)
                c1.metric("결합 영향도", f"{icon} {status}")
                c2.metric("포켓과의 거리", f"{dist} residues")
                pdb_path = get_pdb_file('3WZE')
                if pdb_path:
                    st_molstar(pdb_path, key='her2_viewer', height=400)
else:
    st.error("GDC API 연결 실패")
