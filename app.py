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
    file_path = 'data/clinical.tsv'
    if not os.path.exists(file_path): return None
    
    try:
        df_cli = pd.read_csv(file_path, sep='\t')
        
        # [진단용] 실제 컬럼명과 상위 3개 데이터를 화면에 출력
        st.write("🔍 **임상 데이터 실제 컬럼 목록:**", df_cli.columns.tolist())
        st.write("📊 **데이터 상단 샘플:**", df_cli.head(3))
        
        # HER2 관련 단어가 포함된 모든 컬럼 찾기
        ihc_col = next((c for c in df_cli.columns if 'her2' in c.lower() and 'ihc' in c.lower()), None)
        fish_col = next((c for c in df_cli.columns if 'her2' in c.lower() and 'fish' in c.lower()), None)
        
        # 만약 위 조건으로 못 찾으면 'her2'가 들어간 모든 컬럼이라도 확보
        if not ihc_col:
            ihc_col = next((c for c in df_cli.columns if 'her2' in c.lower()), None)

        def check_her2_low(row):
            val = str(row.get(ihc_col, '')).strip().upper()
            # 1+, 2+, Positive, Low 등 파일에 적힌 실제 값을 확인해야 합니다.
            if val in ['1+', 'IHC 1+', '1']: return True
            # FISH 데이터가 없는 경우를 대비해 IHC 2+만으로도 일단 True로 잡고 테스트
            if val in ['2+', 'IHC 2+', '2']: return True 
            return False

        df_cli['is_her2_low'] = df_cli.apply(check_her2_low, axis=1)
        return df_cli
    except Exception as e:
        st.error(f"진단 중 오류 발생: {e}")
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
    
    # 데이터를 더 많이 가져오기 위해 fields를 확장하고 expand 옵션을 활용합니다.
    params = {
        "filters": json.dumps(filters),
        "fields": "consequence.transcript.aa_change,occurrence.case.submitter_id",
        "format": "JSON",
        "size": "2000" # 2,000개로 대폭 확장
    }
    
    try:
        r = requests.get(ssm_url, params=params)
        res_json = r.json()
        hits = res_json['data']['hits']
        
        data = []
        for h in hits:
            aa = h.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            # 한 변이에 여러 환자(occurrence)가 있을 수 있으므로 모두 추출
            occurrences = h.get('occurrence', [])
            for occ in occurrences:
                case_id = occ.get('case', {}).get('submitter_id')
                data.append({"Case_ID": case_id, "AA_Change": aa})
        
        df = pd.DataFrame(data)
        # 중복 제거 전 데이터 개수 확인용 로그 (Streamlit 콘솔에 찍힘)
        print(f"Total rows fetched: {len(df)}") 
        return df
    except Exception as e:
        st.error(f"데이터 증폭 중 오류: {e}")
        return None

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
# 4. 데이터 통합 및 필터링 적용 부분 수정
if df_mut is not None:
    df_display = df_mut.copy()
    
    if her2_low_only:
        if df_clinical is not None:
            # ID 컬럼 확인 (cases.case_id를 우선적으로 확인)
            id_col = 'cases.case_id' if 'cases.case_id' in df_clinical.columns else \
                     next((col for col in ['case_submitter_id', 'case_id'] if col in df_clinical.columns), None)

            if id_col:
                # [진단 코드 추가] 실제 HER2-Low로 분류된 환자가 있는지 확인
                low_patients = df_clinical[df_clinical['is_her2_low'] == True]
                low_ids = low_patients[id_col].unique()
                
                # 사이드바에 진단 정보 출력
                st.sidebar.info(f"임상 데이터 내 HER2-Low 환자수: {len(low_ids)}명")
                if len(low_ids) > 0:
                    st.sidebar.write("임상 ID 샘플:", list(low_ids)[:3])
                    st.sidebar.write("변이 ID 샘플:", df_mut['Case_ID'].unique()[:3].tolist())

                df_display = df_mut[df_mut['Case_ID'].isin(low_ids)]
                
                if len(df_display) > 0:
                    st.sidebar.success(f"매칭 성공: {len(df_display)}명의 변이 데이터 표시")
                else:
                    st.sidebar.warning("ID 형식이 일치하지 않아 매칭된 데이터가 없습니다.")
            else:
                st.sidebar.error("ID 컬럼을 찾을 수 없습니다.")
        else:
            st.sidebar.error("clinical.tsv 파일을 로드하지 못했습니다.")

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
st.write(df_clinical.columns.tolist())
