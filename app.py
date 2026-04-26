import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정 및 함수 정의
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# [신규] 임상 데이터 로드 및 HER2-Low 분류 함수
@st.cache_data
def load_clinical_data():
    """data/clinical.tsv 파일을 읽어 HER2-Low 환자군을 분류합니다."""
    file_path = 'data/clinical.tsv'
    if not os.path.exists(file_path):
        return None
    
    try:
        # TSV 파일 읽기
        df_cli = pd.read_csv(file_path, sep='\t')
        
        # HER2-Low 정의 로직: IHC 1+ 또는 (IHC 2+ 이면서 FISH Negative)
        def check_her2_low(row):
            # 컬럼명은 GDC 데이터 표준에 따라 'her2_status_by_ihc', 'her2_fish_status' 가정
            # 데이터마다 다를 수 있으므로 .get()으로 안전하게 접근
            ihc = str(row.get('her2_status_by_ihc', '')).strip().upper()
            fish = str(row.get('her2_fish_status', '')).strip().upper()
            
            if ihc == '1+':
                return True
            if ihc == '2+' and (fish == 'NEGATIVE' or fish == 'NON-AMPLIFIED'):
                return True
            return False

        df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
        return df_cli
    except Exception as e:
        st.error(f"임상 데이터 파싱 오류: {e}")
        return None

def estimate_binding_energy(res_num_str, drug_pocket_center=755):
    """변이 지점과 약물 결합 포켓 중심 사이의 거리를 계산"""
    try:
        res_int = int("".join(filter(str.isdigit, str(res_num_str))))
        distance = abs(res_int - drug_pocket_center)
        if distance < 15: return "Critical", "🔴", distance
        elif distance < 30: return "Moderate", "🟡", distance
        else: return "Low", "🟢", distance
    except:
        return "Unknown", "⚪", 0

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
    params = {"filters": json.dumps(filters), "fields": "consequence.transcript.aa_change,occurrence.case.submitter_id", "format": "JSON", "size": "200"}
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

# 2. 메인 화면 구성
st.title("🧬 HER2(ERBB2) Analysis Platform")

# 사이드바 필터 설정
st.sidebar.header("🔍 Filter Settings")
her2_low_only = st.sidebar.checkbox("HER2-Low 환자군만 보기")

with st.spinner('데이터를 분석 중입니다...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 3. 데이터 통합 및 필터링 적용
# 3. 데이터 통합 및 필터링 적용
if df_mut is not None:
    if her2_low_only:
        if df_clinical is not None:
            # [수정] ID 컬럼명을 유연하게 찾기 (GDC 표준 후보군들)
            id_candidates = ['case_submitter_id', 'case_id', 'entity_submitter_id', 'submitter_id']
            id_col = next((col for col in id_candidates if col in df_clinical.columns), None)

            if id_col:
                # HER2-Low 환자의 ID 추출
                low_ids = df_clinical[df_clinical['is_her2_low'] == True][id_col].unique()
                df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
                st.sidebar.success(f"HER2-Low 환자 {len(df_display)}명의 데이터 표시 중")
            else:
                # ID 컬럼을 아예 못 찾은 경우 진단 정보 출력
                st.sidebar.error("ID 컬럼을 찾을 수 없습니다.")
                st.sidebar.write("파일 컬럼 목록:", df_clinical.columns.tolist()[:5]) # 상위 5개만 출력
                df_display = df_mut
        else:
            st.sidebar.error("clinical.tsv 파일을 로드하지 못했습니다.")
            df_display = df_mut
    else:
        df_display = df_mut

    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            top_mats = df_display['AA_Change'].value_counts().reset_index()
            top_mats.columns = ['Mutation', 'Count']
            st.dataframe(top_mats.head(10), use_container_width=True)
        else:
            st.warning("선택한 조건에 해당하는 환자 데이터가 없습니다.")
            top_mats = pd.DataFrame(columns=['Mutation', 'Count'])

    with col2:
        st.subheader("🔬 3D Structure & Energy Analysis")
        
        # 데이터가 있을 때만 분석 진행
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
            st.info("분석할 변이 데이터가 없습니다.")
else:
    st.error("GDC API 연결 실패")
