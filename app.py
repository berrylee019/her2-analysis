import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정 및 함수 정의 (가장 먼저 정의되어야 함)
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

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
    params = {"filters": json.dumps(filters), "fields": "genomic_dna_change,consequence.transcript.aa_change,occurrence.case.submitter_id", "format": "JSON", "size": "100"}
    try:
        r = requests.get(ssm_url, params=params)
        hits = r.json()['data']['hits']
        data = []
        for h in hits:
            aa = h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            data.append({"Case_ID": h.get('occurrence', [{}])[0].get('case', {}).get('submitter_id'), "AA_Change": aa})
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

# [다음 단계 예고] HER2-Low 필터 인터페이스 (기능은 임상 데이터 연결 후 활성화)
st.sidebar.header("Filter Settings")
her2_low_only = st.sidebar.checkbox("HER2-Low 환자군만 보기 (Beta)")
if her2_low_only:
    st.sidebar.info("임상 데이터(IHC/FISH) 연동 작업 중입니다.")

with st.spinner('데이터를 분석 중입니다...'):
    df = get_her2_mutations()

if df is not None:
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.subheader("📊 Mutation Frequency")
        top_mats = df['AA_Change'].value_counts().reset_index()
        top_mats.columns = ['Mutation', 'Count']
        st.dataframe(top_mats.head(10), width='stretch')

    with col2:
        st.subheader("🔬 3D Structure & Energy Analysis")
        selected_mut = st.selectbox("분석할 변이를 선택하세요:", [m for m in top_mats['Mutation'].unique() if m != 'N/A'])
        
        if selected_mut:
            res_num = "".join(filter(str.isdigit, str(selected_mut)))
            status, icon, dist = estimate_binding_energy(res_num)
            
            # 상단 전광판
            c1, c2 = st.columns(2)
            c1.metric("결합 영향도", f"{icon} {status}")
            c2.metric("포켓과의 거리", f"{dist} residues")
            
            # [강력 가이드 로직] 755번일 때만 나타나는 비밀 병기
            if res_num == "755":
                st.error("🚨 [TARGET LOCATED] 755번 변이 지점 정밀 분석 중")
                
                # 3D 뷰어와 가이드 이미지를 나란히 배치
                v_col1, v_col2 = st.columns([2, 1])
                
                with v_col1:
                    pdb_path = get_pdb_file('3WZE')
                    if pdb_path:
                        st_molstar(pdb_path, key='her2_viewer', height=400)
                
                with v_col2:
                    st.write("📍 **위치 가이드**")
                    # 형광색으로 표시된 참고 이미지를 띄워 형님의 눈을 도와드립니다.
                    # (이미지가 없다면 아래 안내 텍스트가 형광색 효과를 대신합니다)
                    st.markdown("""
                    <div style="background-color: #ccff00; padding: 10px; border-radius: 5px; color: black; font-weight: bold; text-align: center;">
                    💡 3D 화면을 돌려<br>가장 깊은 구멍 안쪽<br>바닥을 보세요!
                    </div>
                    """, unsafe_allow_html=True)
                    st.info("755번은 약물이 앉는 '방석'의 정중앙입니다.")
            
            else:
                # 일반 변이일 때는 기존대로 3D만 출력
                pdb_path = get_pdb_file('3WZE')
                if pdb_path:
                    st_molstar(pdb_path, key='her2_viewer', height=400)
