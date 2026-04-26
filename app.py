import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# app.py에 추가할 에너지 추정 로직 (예시)
def estimate_binding_energy(res_num, drug_pocket_center=755):
    """
    변이 지점과 약물 결합 포켓 중심 사이의 거리를 기반으로 
    결합 에너지의 변화 가능성을 추정합니다.
    """
    distance = abs(int(res_num) - drug_pocket_center)
    
    if distance < 5:
        return "Critical (결합력에 직접적 영향 가능성 높음)"
    elif distance < 15:
        return "Moderate (간접적인 구조적 변화 가능성)"
    else:
        return "Low (결합 부위와 거리가 멂)"

# UI 반영
st.subheader("⚡ Drug Binding Energy Estimation")
energy_impact = estimate_binding_energy(res_num)
st.metric(label="변이의 약물 결합 영향도", value=energy_impact)

# 1. 페이지 설정 및 스타일 최적화
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

st.title("🧬 HER2(ERBB2) Analysis Platform")
st.markdown("TCGA-BRCA 데이터를 활용한 HER2 변이 분석 및 3D 매핑 결과입니다.")

# 2. PDB 파일 다운로드 함수 (파일 없음 에러 방지)
def get_pdb_file(pdb_id):
    file_path = f"{pdb_id}.pdb"
    if not os.path.exists(file_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, "w") as f:
                f.write(response.text)
        else:
            st.error(f"PDB 파일을 다운로드할 수 없습니다: {pdb_id}")
            return None
    return file_path

# 3. 데이터 로드 및 분석 함수
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
        "fields": "genomic_dna_change,mutation_subtype,consequence.transcript.aa_change,occurrence.case.submitter_id",
        "format": "JSON",
        "size": "100"
    }
    try:
        r = requests.get(ssm_url, params=params)
        r.raise_for_status()
        hits = r.json()['data']['hits']
        data = []
        for h in hits:
            consq = h.get('consequence', [{}])[0].get('transcript', {})
            data.append({
                "Case_ID": h.get('occurrence', [{}])[0].get('case', {}).get('submitter_id'),
                "AA_Change": consq.get('aa_change', 'N/A'),
                "Type": h.get('mutation_subtype')
            })
        return pd.DataFrame(data)
    except:
        return None

# 4. 메인 실행 로직
with st.spinner('데이터를 분석 중입니다...'):
    df = get_her2_mutations()

if df is not None:
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Mutation Frequency")
        top_mats = df['AA_Change'].value_counts().reset_index()
        top_mats.columns = ['Mutation', 'Count']
        # 경고 해결: width='stretch' 적용
        st.dataframe(top_mats.head(10), width='stretch')

    with col2:
        st.subheader("🔬 3D Structure Mapping")
        selected_mut = st.selectbox("변이를 선택하세요:", top_mats['Mutation'].unique())
        
        # PDB ID '3WZE' (HER2 Kinase domain)
        pdb_path = get_pdb_file('3WZE')
        
        if pdb_path:
            st_molstar(pdb_path, key='her2_viewer', height=400)
            st.caption(f"Target: {selected_mut} on HER2 Structure (PDB: 3WZE)")
else:
    st.error("GDC API 연결에 실패했습니다.")
