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
        # [해결포인트] 사용자가 선택한 후에만 res_num이 정의되도록 로직 구성
        selected_mut = st.selectbox("분석할 변이를 선택하세요:", [m for m in top_mats['Mutation'].unique() if m != 'N/A'])
        
        if selected_mut:
            # 1. 여기서 res_num을 추출합니다.
            res_num = "".join(filter(str.isdigit, str(selected_mut)))
            
            # 2. 추출된 번호로 에너지 계산 함수 호출
            status, icon, dist = estimate_binding_energy(res_num)
            
            # 3. 결과 출력
            c1, c2 = st.columns(2)
            c1.metric("결합 영향도", f"{icon} {status}")
            c2.metric("포켓과의 거리", f"{dist} residues")

            # [강조 로직] 755번인 경우 특별 안내 추가
            if res_num == "755":
                st.warning("🚨 [급소 포착] 755번은 약물 결합 동굴의 핵심 벽면입니다!")
                st.markdown("---")
                st.info("💡 **3D 뷰어 조작 팁:**\n1. 마우스 휠로 **동굴 안쪽**까지 줌인하세요.\n2. 755번은 **동굴 바닥의 오른쪽 벽면**에 위치합니다.")
                
            pdb_path = get_pdb_file('3WZE')
            if pdb_path:
                st_molstar(pdb_path, key='her2_viewer', height=400)
                st.caption(f"📍 분석 지점: {selected_mut} (활성 부위 755번 기준)")
                
                # 뷰어 하단에 형광색 가이드바 표시
                if res_num == "755":
                    st.write("🟢 **형광색 가이드:** 뷰어 안쪽의 가장 깊은 골짜기를 확인하세요.")
else:
    st.error("GDC API 연결 실패")
