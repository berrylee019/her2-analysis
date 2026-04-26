import requests
import pandas as pd
import json
import streamlit as st
from streamlit_molstar import st_molstar

# app.py 상단에 적용
@st.cache_data
def load_gdc_data(query):
    # API 호출 로직
    return data

@st.cache_resource
def load_3d_model(pdb_code):
    # 단백질 모델 로딩 로직 (리소스 소모가 큼)
    return model
    
def get_her2_mutations():
    # 1. API 엔드포인트 설정 (SSMs: Simple Somatic Mutations)
    ssm_methods_url = "https://api.gdc.cancer.gov/ssms"

    # 2. 필터 설정: TCGA-BRCA 프로젝트이면서 유전자 심볼이 ERBB2인 것만
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

    # 3. 가져올 필드 정의 (변이 종류, 아미노산 변화, 영향력 등)
    fields = [
        "genomic_dna_change",
        "mutation_subtype",
        "consequence.transcript.aa_change", # 아미노산 변화 (예: L755S)
        "consequence.transcript.consequence_type",
        "occurrence.case.submitter_id"
    ]
    fields = ",".join(fields)

    # 4. 파라미터 구성
    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "JSON",
        "size": "100"  # 상위 100개 데이터
    }

    # 5. API 호출
    response = requests.get(ssm_methods_url, params=params)
    
    if response.status_code == 200:
        data = response.json()['data']['hits']
        
        # 데이터 정제 (리스트 형태로 변환)
        refined_data = []
        for hit in data:
            aa_change = hit.get('consequence', [{}])[0].get('transcript', {}).get('aa_change', 'N/A')
            refined_data.append({
                "Case_ID": hit.get('occurrence', [{}])[0].get('case', {}).get('submitter_id'),
                "DNA_Change": hit.get('genomic_dna_change'),
                "AA_Change": aa_change,
                "Type": hit.get('mutation_subtype')
            })
        
        return pd.DataFrame(refined_data)
    else:
        print(f"Error: {response.status_code}")
        return None

# 데이터 호출 및 확인
df_her2 = get_her2_mutations()

if df_her2 is not None:
    # 중복 제거 및 빈도순 정렬
    top_mutations = df_her2['AA_Change'].value_counts().reset_index()
    top_mutations.columns = ['Amino_Acid_Change', 'Frequency']
    print("--- HER2(ERBB2) Top Mutations in TCGA-BRCA ---")
    print(top_mutations.head(10))
