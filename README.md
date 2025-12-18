# Quadratic Placer with Recursive Bipartitioning

이 프로그램은 회로 설계(VLSI Design) 자동화 과정 중 **Placement(배치)** 단계를 수행합니다. Quadratic Placer를 사용하여 게이트 간의 연결 길이를 최소화하는 최적 위치를 찾고, 이진 분할(Bipartitioning)을 통해 세부 배치를 반복적으로 수행합니다.

## 핵심 알고리즘 로직

1. **Quadratic Placement (QP):**

   * 게이트 간의 연결을 스프링 시스템(**$\sum w_{ij}(x_i - x_j)^2 + \sum w_{ij}(y_i - y_j)^2$**)으로 모델링합니다. 이를 편미분하여 $Ax=b$ 형태의 방정식을 $x$좌표에 대해 하나, $y$좌표에 대해 하나 얻습니다. $A$ 행렬은 각 gate와 pad의 연결 상태를 이용해 계산할 수 있으며, 이 행렬은 sparse하며 positive semi-definite임이 알려져 있습니다.
     * $A$ 행렬의 대각 성분은 해당하는 gate의 모든 weight의 합입니다.
     * $A$ 행렬의 대각 성분이 아닌 성분은 두 개의 gate 사이의 weight의 합의 $(-1)$배입니다.
     * 이때 sparse matrix solver를 이용하기 때문에 행렬의 indexing에 주의해야 합니다. index가 0부터 시작하여야 하므로 각 QP 단계에 대한 `start`와 `end`, `num_gates` 배열을 만들어 관리하였습니다. 또한 net에 저장된 gate id와 index를 연결하기 위해 배열(`inverse_gates`)을 만들어 관리하였습니다.
   * 희소 행렬(Sparse Matrix) solver를 사용하여 에너지가 최소화되는 평형 상태의 좌표를 계산합니다.
2. **Weighting Models:**

   * **Default:** 기본 가중치 적용.
   * **Clique Model:** 넷(Net)의 크기에 반비례하는 가중치를 적용하여 다중 연결망의 최적화를 돕습니다.
3. **Fixed Boundary Constraints:**

   * 분할 시 현재 영역 밖에 위치하게 된 게이트나 패드(Pad)들은 해당 영역의 경계선에 고정된 **Pseudo Pad**로 간주되어 힘의 평형 계산에 포함됩니다.
4. **Recursive Bipartitioning:** `exponential_iter+1` 횟수만큼 분할을 반복합니다. 즉, `exponential_iter`가 1일 때 QP1까지 진행하며, 2일때 QP2, 3까지 진행하고, 3일때 QP4, 5, 6, 7까지 진행합니다.

   * 홀수 번째 반복(**$ei=1, 3, \dots$**): solving이 완료된 후 X-좌표 기준으로 정렬 후 좌/우 분할.
   * 짝수 번째 반복(**$ei=2, 4, \dots$**): Y-좌표 기준으로 정렬 후 상/하 분할.
   * 이때 각 QP 단계가 완료되면 다음 QP(각각 2QP와 2QP+1)을 준비하기 위해 해당 구간의 gate들을 정렬합니다. 또한 다음 QP의 물리적 border와 index border를 계산합니다.

---

## 입출력 형식

### 1. 입력 파일 (실행 시 위치 확인)

입력 데이터는 다음의 순서를 따릅니다.

```
<num_gates> <num_nets>
<gate_id> <num_connected_nets> <net_id_1> <net_id_2> ... (num_gates 만큼 반복)
<num_pads>
<pad_id> <connected_net_id> <x_coord> <y_coord> (num_pads 만큼 반복)
```

### 2. 출력 파일 (`final_output.txt`)

배치가 완료된 후 각 게이트의 최종 좌표를 출력합니다.

```
<gate_id> <final_x_coord> <final_y_coord>
```

---

## 실행 방법

### 컴파일

`solver.cpp` 파일이 아래와 같이 지정된 경로에 있어야 합니다. `solver.h` 파일은 `solver.cpp` 파일과 동일한 경로에 있어야 합니다. 

```
g++ -o QuadraticPlacer QuadraticPlacer.cpp ./solvers/cpp/solver.cpp 
```

### 실행

```
# 기본 모드
./QuadraticPlacer.exe <input_file>

# Clique Weight 모드
./QuadraticPlacer.exe <input_file> --clique
```

실행 후 콘솔에 `Enter number of exponential iterations`가 표시되면 원하는 분할 깊이(예: 1, 2, 3...)를 입력합니다.
