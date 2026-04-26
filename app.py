import streamlit as st
import pandas as pd
import requests
import json
import os
from streamlit_molstar import st_molstar

# 1. 페이지 설정
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬", layout="wide")

# 2. 데이터 로드 및 분석 함수
@st.cache_data
def load_clinical_data():
    """GDC CSV의 복잡한 컬럼 구조를 파싱하고 데이터 부재 시 진단 정보를 제공합니다."""
    paths = [
        'data/2026-04-26T13-44_export.csv',
        '2026-04-26T13-44_export.csv',
        'data/clinical.tsv',
        'clinical.tsv'
    ]
    
    target_path = next((p for p in paths if os.path.exists(p)), None)
    if not target_path:
        return None
    
    try:
        # [해결] engine='python' 사용 시 low_memory=False 옵션은 제거해야 에러가 발생하지 않습니다.
        df_cli = pd.read_csv(target_path, sep=None, engine='python')
        
        cols = df_cli.columns.tolist()
        
        # IHC 컬럼 탐색: 패턴 매칭 강화
        ihc_candidates = [c for c in cols if 'her2' in c.lower() and 'ihc' in c.lower()]
        if not ihc_candidates:
            ihc_candidates = [c for c in cols if 'her2_status' in c.lower() and 'fish' not in c.lower()]
        
        ihc_col = ihc_candidates[0] if ihc_candidates else None
        
        # FISH 컬럼 탐색
        fish_candidates = [c for c in cols if 'her2' in c.lower() and 'fish' in c.lower()]
        fish_col = fish_candidates[0] if fish_candidates else None
        
        # ID 컬럼 탐색 (cases.submitter_id 확인용)
        id_col = next((c for c in cols if 'submitter_id' in c.lower()), None)

        # 진단 모드: 필수 데이터(IHC)가 없는 경우
        if not ihc_col:
            return {"status": "MISSING_COLUMNS", "columns": cols}

        def check_her2_low(row):
            ihc_val = str(row.get(ihc_col, '')).strip().upper()
            fish_val = str(row.get(fish_col, '')).strip().upper() if fish_col else ""
            
            # HER2-Low 판정 로직 (숫자 및 텍스트 대응)
            if any(x in ihc_val for x in ['1+', '1.0', '1']):
                return True
            if any(x in ihc_val for x in ['2+', '2.0', '2']):
                neg_terms = ['NEG', 'NON', 'NOT', '0']
                if any(t in fish_val for t in neg_terms): return True
                if fish_val in ['', 'NAN', '--', "'--", 'UNKNOWN', 'PND']: return True
            return False

        if id_col:
            df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
            df_cli['Match_ID'] = df_cli[id_col].astype(str).str.strip()
            return {"status": "SUCCESS", "data": df_cli}
        
        return None

    except Exception as e:
        st.sidebar.error(f"⚠️ 파일 분석 오류: {e}")
        return None

# --- GDC 변이 데이터 호출 ---
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
        data = [{"Case_ID": str(occ.get('case', {}).get('submitter_id')).strip(), 
                 "AA_Change": h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')}
                for h in hits for occ in h.get('occurrence', [])]
        return pd.DataFrame(data)
    except: return None

def get_pdb_file(pdb_id):
    file_path = f"{pdb_id}.pdb"
    if not os.path.exists(file_path):
        try:
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            r = requests.get(url)
            with open(file_path, "w") as f: f.write(r.text)
        except: return None
    return file_path

# 3. 메인 UI
st.title("🧬 HER2(ERBB2) Analysis Platform")

st.sidebar.header("🔍 Filter Settings")
low_only = st.sidebar.checkbox("HER2-Low 환자군 분석")

with st.spinner('데이터 분석 및 매칭 중...'):
    df_mut = get_her2_mutations()
    clinical_result = load_clinical_data()

# 4. 결과 출력 및 진단 가이드
if df_mut is not None:
    df_display = df_mut.copy()
    
    if low_only:
        if clinical_result and clinical_result.get("status") == "SUCCESS":
            df_clinical = clinical_result["data"]
            low_ids = df_clinical[df_clinical['is_her2_low'] == True]['Match_ID'].unique()
            st.sidebar.success(f"✅ HER2-Low: {len(low_ids)}명 식별")
            df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
            st.sidebar.info(f"🎯 매칭된 변이: {len(df_display)}건")
        
        elif clinical_result and clinical_result.get("status") == "MISSING_COLUMNS":
            st.error("⚠️ 업로드된 파일에 **HER2(IHC) 정보**가 포함되어 있지 않습니다.")
            with st.expander("현재 파일에서 확인된 컬럼 목록 (GDC 재다운로드 시 참고)"):
                st.write(clinical_result["columns"])
            st.info("💡 **해결책**: GDC Data Portal에서 Clinical 데이터를 받을 때, **'her2_status_by_ihc'** 항목이 포함되도록 Custom Download를 진행해 주세요.")
            st.stop()
        else:
            st.sidebar.error("임상 데이터를 로드할 수 없습니다.")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("📊 Mutation Frequency")
        if not df_display.empty:
            counts = df_display['AA_Change'].value_counts().reset_index()
            counts.columns = ['Mutation', 'Count']
            st.dataframe(counts.head(10), use_container_width=True)
        else:
            st.warning("분석할 데이터가 없습니다.")

    with col2:
        st.subheader("🔬 3D Structure Analysis")
        if not df_display.empty:
            mats = [m for m in df_display['AA_Change'].unique() if m != 'N/A']
            selected = st.selectbox("변이 선택:", mats)
            if selected:
                res_num_str = "".join(filter(str.isdigit, str(selected)))
                if res_num_str:
                    dist = abs(int(res_num_str) - 755)
                    st.metric("Dist from Pocket (755)", f"{dist} residues")
                
                pdb_path = get_pdb_file('3WZE')
                if pdb_path:
                    st_molstar(pdb_path, height=400)
                else:
                    st.info("PDB 구조 파일을 불러올 수 없습니다.")
else:
    st.error("GDC API 연결 실패")
