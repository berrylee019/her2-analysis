import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 함수 (실패 원천 차단 버전)
@st.cache_data
def load_clinical_data():
    """모든 경로와 파일 형식을 뒤져서 임상 데이터를 반드시 찾아냅니다."""
    # 1. 탐색할 파일 후보들 (형님이 말씀하신 tsv와 업로드된 csv 모두 포함)
    possible_files = [
        'data/clinical.tsv',
        'data/2026-04-26T13-44_export.csv',
        'clinical.tsv',
        '2026-04-26T13-44_export.csv'
    ]
    
    target_path = None
    for p in possible_files:
        if os.path.exists(p):
            target_path = p
            break
            
    if not target_path:
        # 최후의 수단: data 폴더 내의 모든 csv/tsv 검색
        if os.path.exists('data'):
            all_files = [os.path.join('data', f) for f in os.listdir('data') if f.endswith(('.csv', '.tsv'))]
            if all_files: target_path = all_files[0]

    if not target_path: return None

    try:
        # 구분자 자동 판단 (tsv는 \t, csv는 ,)
        sep = '\t' if target_path.endswith('.tsv') else ','
        df_cli = pd.read_csv(target_path, sep=sep)
        
        # [핵심] GDC 특유의 복잡한 컬럼명에서 키워드로 낚아채기
        # 'her2'와 'ihc'가 동시에 들어간 컬럼을 찾습니다.
        ihc_col = next((c for c in df_cli.columns if 'her2' in c.lower() and 'ihc' in c.lower()), None)
        fish_col = next((c for c in df_cli.columns if 'her2' in c.lower() and 'fish' in c.lower()), None)
        # ID 컬럼: 'submitter_id'가 포함된 컬럼 탐색
        id_col = next((c for c in df_cli.columns if 'submitter_id' in c.lower()), None)

        def check_her2_low(row):
            ihc = str(row.get(ihc_col, '')).strip().upper() if ihc_col else ""
            fish = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            # HER2-Low 판정 로직
            if '1+' in ihc: return True
            if '2+' in ihc:
                if any(x in fish for x in ['NEG', 'NON', 'NOT']): return True
                if fish in ['', 'NAN', '--', "'--"]: return True # 데이터가 없는 경우도 포함
            return False

        if ihc_col and id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return df_cli
        return None
    except:
        return None

# ... (중략: get_her2_mutations, estimate_binding_energy 함수는 기존과 동일) ...

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

# 3. 메인 화면
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
her2_low_only = st.sidebar.checkbox("HER2-Low 환자군만 보기")

with st.spinner('데이터를 분석 중입니다...'):
    df_mut = get_her2_mutations()
    df_clinical = load_clinical_data()

# 4. 데이터 통합
if df_mut is not None:
    df_display = df_mut.copy()
    
    if her2_low_only:
        if df_clinical is not None:
            low_patients = df_clinical[df_clinical['is_her2_low'] == True]
            low_ids = low_patients['Match_ID'].unique()
            # [디버깅] 화면에 환자수 표시
            st.sidebar.info(f"HER2-Low 분류 환자: {len(low_ids)}명")
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
        else:
            st.sidebar.error("데이터 로드 실패: 컬럼명을 확인하세요.")

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
            selected_mut = st.selectbox("변이 선택:", mats)
            if selected_mut:
                res_num = "".join(filter(str.isdigit, str(selected_mut)))
                status, icon, dist = estimate_binding_energy(res_num)
                st.metric("결합 영향도", f"{icon} {status}")
                pdb_path = get_pdb_file('3WZE')
                if pdb_path: st_molstar(pdb_path, height=400)
else:
    st.error("GDC API 연결 실패")
