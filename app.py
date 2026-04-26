import requests
import pandas as pd
import json
import streamlit as st
from streamlit_molstar import st_molstar

# 1. 페이지 설정 (가장 상단에 위치해야 함)
st.set_page_config(page_title="HER2 Analysis Platform", page_icon="🧬")

st.title("🧬 HER2(ERBB2) Analysis Platform")
st.markdown("TCGA-BRCA 데이터를 활용한 HER2 변이 분석 결과입니다.")

# 2. 캐싱된 데이터 호출 함수
@st.cache_data
def get_her2_mutations():
    # API 엔드포인트 설정 (SSMs: Simple Somatic Mutations)
    ssm_methods_url = "https://api.gdc.cancer.gov/ssms"

    # 필터 설정: TCGA-BRCA 프로젝트이면서 유전자 심볼이 ERBB2인 것만
    filters = {
        "op": "and",
        "content": [
            {
                "op": "in",
                "content": {
                    "field": "occurrence.case.project.project_id",
                    "value": ["TCGA-BRCA"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "genes.symbol",
                    "value": ["ERBB2"]
                }
            }
        ]
    }

    fields = [
        "genomic_dna_change",
        "mutation_subtype",
        "consequence.transcript.aa_change",
        "consequence.transcript.consequence_type",
        "occurrence.case.submitter_id"
    ]
    
    params = {
        "filters": json.dumps(filters),
        "fields": ",".join(fields),
        "format": "JSON",
        "size": "100"
    }

    try:
        response = requests.get(ssm_methods_url, params=params)
        response.raise_for_status()
        data = response.json()['data']['hits']
        
        refined_data = []
        for hit in data:
            # 안전하게 아미노산 변화 데이터 추출
            consequences = hit.get('consequence', [])
            aa_change = 'N/A'
            if consequences:
                aa_change = consequences[0].get('transcript', {}).get('aa_change', 'N/A')
            
            refined_data.append({
                "Case_ID": hit.get('occurrence', [{}])[0].get('case', {}).get('submitter_id'),
                "DNA_Change": hit.get('genomic_dna_change'),
                "AA_Change": aa_change,
                "Type": hit.get('mutation_subtype')
            })
        
        return pd.DataFrame(refined_data)
    except Exception as e:
        st.error(f"데이터를 가져오는 중 오류 발생: {e}")
        return None

# 3. 메인 로직 실행
with st.spinner('GDC API로부터 데이터를 불러오는 중...'):
    df_her2 = get_her2_mutations()

if df_her2 is not None:
    # 빈도 분석
    top_mutations = df_her2['AA_Change'].value_counts().reset_index()
    top_mutations.columns = ['Amino_Acid_Change', 'Frequency']
    
    # 웹 화면 출력
    st.subheader("📊 Top 10 Mutations in HER2 (TCGA-BRCA)")
    st.dataframe(top_mutations.head(10), use_container_width=True)
    
    # 상세 데이터 확인
    with st.expander("전체 원본 데이터 보기"):
        st.write(df_her2)
else:
    st.warning("데이터를 불러올 수 없습니다. API 상태를 확인해주세요.")

# 4. 3D 모델 로딩 (추후 구현을 위한 플레이스홀더)
st.divider()
st.subheader("🔬 3D Structure Analysis")
st.info("선택한 변이에 따른 단백질 구조 변화를 시각화 합니다.")
# ... (기존 변이 분석 코드 아래에 추가)

if df_her2 is not None:
    st.divider()
    st.subheader("🔬 Mutation 3D Mapping")
    
    # 1. 분석된 변이 중 시각화할 항목 선택
    mutation_list = top_mutations['Amino_Acid_Change'].tolist()
    selected_mut = st.selectbox("구조를 확인할 변이를 선택하세요:", [m for m in mutation_list if m != 'N/A'])

    if selected_mut:
        # 2. 변이 번호 추출 (예: L755S -> 755)
        try:
            res_num = "".join(filter(str.isdigit, selected_mut))
            wild_aa = selected_mut[0]
            mut_aa = selected_mut[-1]
            
            st.info(f"선택된 변이: {wild_aa} (Wild-type) ➔ {res_num}번 위치 ➔ {mut_aa} (Mutant)")

            # 3. 3D 시각화 (HER2 Kinase Domain 구조: 3WZE 사용)
            # st_molstar는 PDB ID나 URL을 직접 받을 수 있습니다.
            # '3WZE'는 HER2의 활성 상태 구조입니다.
            
            col1, col2 = st.columns([3, 1])
            
            with col1:
                st_molstar(
                    '3WZE', 
                    key='her2_viewer',
                    height=500,
                )
            
            with col2:
                st.write("**분석 가이드**")
                st.caption(f"""
                - **PDB ID:** 3WZE
                - **Target Residue:** {res_num}
                - 이 부위가 약물(Lapatinib 등) 결합 포켓 근처인지 확인하십시오.
                - {res_num}번 잔기의 부피 변화가 구조 안정성에 미치는 영향을 분석할 수 있습니다.
                """)
                
        except Exception as e:
            st.error(f"구조 매핑 중 오류가 발생했습니다: {e}")

# 3D 구조와 약물 결합 이해를 돕기 위한 이미지
#
