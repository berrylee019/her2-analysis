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
    # 파일 탐색 (형님이 업로드하신 파일명 최우선)
    paths = [
        'data/2026-04-26T13-44_export.csv',
        'data/clinical.tsv',
        '2026-04-26T13-44_export.csv',
        'clinical.tsv'
    ]
    
    target_path = next((p for p in paths if os.path.exists(p)), None)
    if not target_path:
        return None
    
    try:
        sep = '\t' if target_path.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_path, sep=sep)
        
        # [해결책] 패턴 매칭 방식으로 컬럼 찾기
        # GDC CSV는 'diagnoses.0.her2_status_by_ihc' 처럼 중간에 숫자가 들어갑니다.
        cols = df_cli.columns.tolist()
        
        # IHC 컬럼 찾기: 'her2'와 'ihc'가 동시에 포함된 첫 번째 컬럼
        ihc_col = next((c for c in cols if 'her2' in c.lower() and 'ihc' in c.lower()), None)
        # FISH 컬럼 찾기: 'her2'와 'fish'가 동시에 포함된 첫 번째 컬럼
        fish_col = next((c for c in cols if 'her2' in c.lower() and 'fish' in c.lower()), None)
        # ID 컬럼 찾기: 'submitter_id'가 포함된 컬럼
        id_col = next((c for c in cols if 'submitter_id' in c.lower()), 
                      next((c for c in cols if 'cases.case_id' in c.lower()), None))

        # 만약 IHC 컬럼을 못 찾았다면 'her2_status'만 들어간 거라도 확보
        if not ihc_col:
            ihc_col = next((c for c in cols if 'her2_status' in c.lower() and 'fish' not in c.lower()), None)

        def check_her2_low(row):
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 판정 로직
            # 1) IHC 1+ 면 무조건 Low
            if '1+' in ihc or ihc == '1': return True
            # 2) IHC 2+ 일 때는 FISH 음성인 경우만 Low
            if '2+' in ihc or ihc == '2':
                if any(x in fish for x in ['NEG', 'NON', 'NOT']): return True
                if fish in ['', 'NAN', '--', "'--", 'UNKNOWN']: return True
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        else:
            # 실패 시 디버깅을 위해 사이드바에 정보 노출
            st.sidebar.warning(f"매칭 시도 컬럼 - IHC: {ihc_col}, ID: {id_col}")
            return None
    except Exception as e:
        st.sidebar.error(f"파일 처리 중 오류: {e}")
        return None

# ... (중략: get_her2_mutations, estimate_binding_energy, get_pdb_file 함수는 기존과 동일) ...
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

def estimate_binding_energy(res_num_str, drug_pocket_center=755):
    try:
        res_int = int("".join(filter(str.isdigit, str(res_num_str))))
        distance = abs(res_int - drug_pocket_center)
        if distance < 15: return "Critical", "🔴", distance
        elif distance < 30: return "Moderate", "🟡", distance
        else: return "Low", "🟢", distance
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
            
            # GDC 변이 데이터와 매칭
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            if not df_display.empty:
                st.sidebar.success(f"매칭 성공: {len(df_display)}건")
            else:
                st.sidebar.warning("ID 매칭 실패 (ID 형식을 확인하세요)")
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
                res_num = "".join(filter(str.isdigit, str(selected_mut)))
                status, icon, dist = estimate_binding_energy(res_num)
                st.metric("결합 영향도", f"{icon} {status}")
                pdb_path = get_pdb_file('3WZE')
                if pdb_path: st_molstar(pdb_path, height=400)
else:
    st.error("GDC API 연결 실패")
