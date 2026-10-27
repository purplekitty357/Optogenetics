// 광유전학 1형 당뇨병 치료 시뮬레이터 - 계산 모델
// 모든 식은 "광유전학을 활용한 1형 당뇨병 치료 시뮬레이터 - 실제 수식과 임의 변형 수식 정리" 문서의 번호를 따른다.
//   1.1 Bergman minimal model (실제)        2.1 빛 활성화 dA/dt (임의 변형)
//   1.2 식사 후 포도당 흡수 UG(t) (실제)     2.2 빛-인슐린 분비 (임의 변형)
//   1.3 피하주사 2구획 흡수 (실제)           2.3 Hill 함수 분비 포화 (임의 변형)
//                                            2.4 단순 혈당식 (임의 변형)
//                                            2.5 피하주사 1구획 흡수 (임의 변형)
// 시간 단위: 분(min), 혈당: mg/dL, 인슐린: μU/mL
// 주의: 문서에는 파라미터 수치가 없으므로 아래 기본값은 시연용으로 정한 값이다.

const DEFAULT_PARAMS = {
    // 시뮬레이션
    duration: 18 * 60,  // 총 시간 (분)
    dt: 0.5,            // 반복 계산 간격 (분)
    startHour: 6,       // 시작 시각
    G0: 180,            // 초기 혈당

    // 식사 (1.2)
    meals: [            // 시작 후 몇 분에 몇 g 탄수화물
        { t: 60, carbs: 60 },
        { t: 360, carbs: 70 },
        { t: 720, carbs: 80 },
    ],
    AG: 0.8,            // 탄수화물 생체이용률
    tmaxG: 40,          // 포도당 출현률 최대 시간 (분)
    VG: 160,            // 포도당 분포 부피 (dL) - mg/min 을 mg/dL/min 으로 바꾸는 데 사용

    // 혈당식 선택: 'bergman' (1.1) 또는 'simple' (2.4)
    glucoseModel: 'bergman',
    // 1.1 Bergman
    p1: 0.005,          // 인슐린과 무관한 혈당 제거율 (1형 당뇨는 작음)
    Gb: 150,            // 기저 혈당
    p2: 0.025,          // 인슐린 작용 소멸률
    p3: 1.3e-5,         // 인슐린 -> 작용 X 전환 비율
    n: 0.093,           // 인슐린 제거율
    Ib: 0,              // 기저 인슐린 (β세포 파괴 -> 0)
    // 2.4 단순식
    kbasal: 0.0008,     // 인슐린과 무관한 기본 혈당 제거율
    kinsulin: 0.0005,   // 인슐린의 혈당 저하 효과 크기

    // 2.1 빛 활성화
    kon: 0.5,           // 빛 켤 때 활성화 속도
    koff: 0.1,          // 빛 끌 때 비활성화 속도
    lightDelay: 0,      // 식사 후 몇 분 뒤 빛을 켤지
    lightDuration: 120, // 빛을 켜 두는 시간 (분)

    // 2.2 / 2.3 빛-인슐린 분비
    secretionModel: 'glucose', // 'linear': α·A, 'glucose': α·A·G, 'hill': α·A·SG(G)
    uBasal: 0,          // 기저 인슐린 분비
    alphaLinear: 2.5,   // α (α·A 형태)
    alphaGlucose: 0.015,// α (α·A·G 형태)
    alphaHill: 1.0,     // α (α·A·SG(G) 형태)
    Smax: 4,            // 최대 인슐린 분비율
    G50: 150,           // 최대 분비의 절반을 만드는 혈당
    hillN: 4,           // Hill 지수

    // 주사 모델 선택: 'two' (1.3) 또는 'one' (2.5)
    injectionModel: 'two',
    insulinPerCarb: 1 / 12, // 탄수화물 12 g 당 1 U
    injDelay: 0,        // 식사 몇 분 뒤에 주사하는지
    VI: 12,             // 인슐린 분포 부피 (L): 1 mU/L = 1 μU/mL
    tmaxI: 55,          // 1.3 피하 흡수 최대 시간
    ka: 0.02,           // 2.5 피하 흡수율
};

// 1.2 식사 후 포도당 유입률 UG(t) [mg/min], 여러 끼니를 합한다
function mealAppearance(t, p) {
    let total = 0;
    for (const meal of p.meals) {
        const s = t - meal.t;
        if (s <= 0) continue;
        const DG = meal.carbs * 1000; // g -> mg
        total += DG * p.AG * s * Math.exp(-s / p.tmaxG) / (p.tmaxG * p.tmaxG);
    }
    return total;
}

// 빛 자극 L(t): 각 식사 후 lightDelay 분부터 lightDuration 분 동안 켜짐 (1), 그 외 꺼짐 (0)
function lightOn(t, p) {
    for (const meal of p.meals) {
        const start = meal.t + p.lightDelay;
        if (t >= start && t < start + p.lightDuration) return 1;
    }
    return 0;
}

// 2.3 Hill 함수
function hill(G, p) {
    const Gn = Math.pow(Math.max(G, 0), p.hillN);
    return p.Smax * Gn / (Math.pow(p.G50, p.hillN) + Gn);
}

// 2.2 / 2.3 빛 자극에 의한 인슐린 분비 입력 u_light(t)
function lightSecretion(A, G, p) {
    if (p.secretionModel === 'linear') return p.uBasal + p.alphaLinear * A;
    if (p.secretionModel === 'hill') return p.uBasal + p.alphaHill * A * hill(G, p);
    return p.uBasal + p.alphaGlucose * A * G;
}

// 한 가지 치료 방식으로 시간에 따라 반복 계산 (오일러 방법)
// mode: 'none' (무치료), 'light' (광유전학), 'injection' (주사)
function simulate(mode, p) {
    const steps = Math.round(p.duration / p.dt);
    let G = p.G0, X = 0, I = p.Ib, A = 0;
    let S1 = 0, S2 = 0, Ssc = 0;
    const doses = p.meals.map(m => ({ t: m.t + p.injDelay, units: m.carbs * p.insulinPerCarb, given: false }));
    const out = { t: [], G: [], I: [], A: [], L: [], u: [], doses: [] };

    for (let k = 0; k <= steps; k++) {
        const t = k * p.dt;

        // 주사: δ(t - tinj) 는 해당 시점에 피하 구획에 Dinj 를 한 번에 더하는 것으로 처리
        let L = 0, u = 0;
        if (mode === 'injection') {
            for (const d of doses) {
                if (!d.given && t >= d.t) {
                    const mU = d.units * 1000 / p.VI; // U -> μU/mL 단위의 양
                    if (p.injectionModel === 'two') S1 += mU; else Ssc += mU;
                    d.given = true;
                    out.doses.push({ t, units: d.units });
                }
            }
            u = p.injectionModel === 'two' ? S2 / p.tmaxI : p.ka * Ssc;
        } else if (mode === 'light') {
            L = lightOn(t, p);
            u = lightSecretion(A, G, p);
        }

        out.t.push(t); out.G.push(G); out.I.push(I); out.A.push(A); out.L.push(L); out.u.push(u);

        // 미분값 계산
        const Rfood = mealAppearance(t, p) / p.VG; // mg/dL/min
        let dG, dX = 0, dI;
        if (p.glucoseModel === 'simple') {
            // 2.4  dG/dt = Rfood - kbasal·G - kinsulin·I·G
            dG = Rfood - p.kbasal * G - p.kinsulin * I * G;
        } else {
            // 1.1  dG/dt = -p1(G - Gb) - X·G  (+ 식사 유입 Rfood)
            dG = Rfood - p.p1 * (G - p.Gb) - X * G;
            dX = -p.p2 * X + p.p3 * (I - p.Ib);
        }
        if (mode === 'injection' && p.injectionModel === 'one') {
            // 2.5  dI/dt = ka·Ssc - ke·I  (ke = n)
            dI = u - p.n * I;
        } else {
            // 1.1  dI/dt = -n(I - Ib) + u(t)
            dI = -p.n * (I - p.Ib) + u;
        }
        const dA = p.kon * L * (1 - A) - p.koff * A;          // 2.1
        const dS1 = -S1 / p.tmaxI;                            // 1.3 (u 는 주사 시점에 S1 에 더함)
        const dS2 = S1 / p.tmaxI - S2 / p.tmaxI;              // 1.3
        const dSsc = -p.ka * Ssc;                             // 2.5

        G = Math.max(G + dG * p.dt, 0);
        X = X + dX * p.dt;
        I = Math.max(I + dI * p.dt, 0);
        A = Math.min(Math.max(A + dA * p.dt, 0), 1);
        S1 += dS1 * p.dt; S2 += dS2 * p.dt; Ssc += dSsc * p.dt;
    }
    return out;
}

// 결과 요약 지표
function summarize(r) {
    const n = r.G.length;
    let sum = 0, max = -Infinity, min = Infinity, inRange = 0, low = 0, high = 0;
    for (const g of r.G) {
        sum += g; max = Math.max(max, g); min = Math.min(min, g);
        if (g < 70) low++; else if (g > 180) high++; else inRange++;
    }
    return { mean: sum / n, max, min, tir: 100 * inRange / n, tbr: 100 * low / n, tar: 100 * high / n };
}

if (typeof module !== 'undefined') module.exports = { DEFAULT_PARAMS, simulate, summarize, mealAppearance, hill };
