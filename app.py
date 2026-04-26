import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. [중요] 페이지 설정은 반드시 모든 'st' 코드 중 가장 상단에 위치해야 합니다.
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수 정의
@st.cache_data
def load_clinical_data():
    """data/clinical.tsv 파일을 읽어 HER2-Low 환자군을 분류합니다."""
    file_path = 'data/clinical.tsv'
    if not os.path.exists(file_path):
        return None
    try:
        df_cli = pd.read_csv(file_path, sep='\t')
        # GDC의 'cases.' 접두사가 붙은 컬럼명 대응
        ihc_col = next((c for c in df_cli.columns if 'her2_status_by_ihc' in c), None)
        fish_col = next((c for c in df_cli.columns if 'her2_fish_status' in c), None)

        def check_her2_low(row):
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            if ihc == '1+': return True
            if ihc == '2+' and (fish == 'NEGATIVE' or fish == 'NON-AMPLIFIED'): return True
            return False

        df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
        return df_cli
    except Exception as e:
        st.error(f"임상 데이터 파싱 오류: {e}")
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
    params = {"filters": json.dumps(filters), "fields": "consequence.transcript.aa_change,occurrence.case.submitter_id", "format": "JSON", "size": "500"}
    try:
        r = requests.get(ssm_url, params=params)
        hits = r.json()['data']['hits']
        data = []
        for h in hits:
            aa = h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            case_id = h.get('occurrence', [{}])[0].get('case', {}).get('submitter_id')
            data.append({"Case_ID": case_id, "AA_Change": aa})
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

# 4. 데이터 통합 및 필터링 적용
if df_mut is not None:
    df_display = df_mut.copy()
    
    if her2_low_only:
        if df_clinical is not None:
            # 형님이 확인해주신 'cases.case_id' 컬럼 사용
            id_col = 'cases.case_id' if 'cases.case_id' in df_clinical.columns else \
                     next((col for col in ['case_submitter_id', 'case_id'] if col in df_clinical.columns), None)

            if id_col:
                low_ids = df_clinical[df_clinical['is_her2_low'] == True][id_col].unique()
                df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
                st.sidebar.success(f"HER2-Low 환자 {len(df_display)}명 필터링 완료")
            else:
                st.sidebar.error("ID 컬럼을 찾을 수 없습니다.")
        else:
            st.sidebar.error("clinical.tsv 파일을 찾을 수 없습니다.")

    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            top_mats = df_display['AA_Change'].value_counts().reset_index()
            top_mats.columns = ['Mutation', 'Count']
            st.dataframe(top_mats.head(10), use_container_width=True)
        else:
            st.warning("데이터가 없습니다.")
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
