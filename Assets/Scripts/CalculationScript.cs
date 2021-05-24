using System;
using System.Collections;
using System.Collections.Generic;
// using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Threading;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.UI;
using Debug = UnityEngine.Debug;
// using Random = UnityEngine.Random;

public class IonsParameters {
    public double q; 
    public double m;
    public Color col;
}

public class CalculationScript: MonoBehaviour
{
    private static CultureInfo Culture = CultureInfo.InvariantCulture;
    
    public bool calculatingSimplePowerLines;
    public Toggle simplePowerLinesToggle;

    public bool calculatingPowerLines;
    public Toggle powerLinesToggle;
    
    public bool calculatingIonTrajectory;
    public Toggle ionTrajectoryToggle;
    
    public bool   mField;
    public Toggle mFieldToggle;
    
    public bool   mPlusEFields;
    public Toggle mPlusEFieldsToggle;
    
    public bool   mPlusDEFields;
    public Toggle mPlusDEFieldsToggle;

    public ParticleSystem fieldPoints;
    public List<ParticleSystem.Particle> particleList = new List<ParticleSystem.Particle>();

    private DVector3 B = new DVector3();
    public DQuaternionList MList = new DQuaternionList();

    public Material matMag;
    public Material matLine1;
    public Material matLine2;
    public Material matLine3;

    public float colorIntensity;

    private double z; // м
    private double R; // м, радиус контура S
    private double i; // А, величина силы тока, протекающего по контуру
    // public double O = 0.0f; // центр контура S
    private double ro; // м, расстояние от точки M до оси Oz
    // public double dl = 0.0f; // элементарный участок контура
    public double doubleAlpha; //град, угол между двумя элементарными участками контура, расположенными симметрично
    private double beta;
    private double r; //м, расстояние от точки M до элементарных участков контура dl
    public double mu0 = 1.256637 * Math.Pow(10, -6); //1,256637⋅10⁻⁶ Н·А⁻², магнитная постоянная

    private double K;
    private double L;
    private double B_ro;
    private double B_z;

    public double maxM;
    public double minM;

    public double meshStep = 0.003; // м, шаг обсчёта
    public InputField ui_meshStep;

    public GameObject lineSource;
    public GameObject fieldSource;
    public GameObject drawZone;

    private GameObject[] lineArray;
    private Vector3[] magneticFieldArray;

    // Задание параметров рассматриваемой области
    public float xFieldSize;
    public float yFieldSize;
    public float zFieldSize;
    public InputField ui_xFieldSize, ui_yFieldSize, ui_zFieldSize;

    public int maxSteps = 200;
    public InputField ui_maxSteps;

    public Transform coilsRoot;
    private List<Coil> coilsP = new List<Coil>();

    public double h; // Масштаб по времени
    public InputField ui_h;
    
    // public string ionName;
    public int ionSteps;
    public InputField ui_ionSteps;

    public Vector3 p0;
    public InputField p0_X, p0_Y, p0_Z;
    public Vector3 p1;
    public InputField p1_X, p1_Y, p1_Z;
    
    public float gammaAngle;
    public float betaAngle;
    float maxIonRange = 0;
    
    // public static double firstIonQ = -1.60217 * Math.Pow(10, -19); //Кл, Electron
    // public static double firstIonM = 9.1 * Math.Pow(10, -31); //кг
    public static double firstIonQ = 1.60217 * Math.Pow(10, -19); //Кл, Proton
    public static double firstIonM = 1.6 * Math.Pow(10, -27); //кг
    public static double secondIonQ = 1.60217 * Math.Pow(10, -19); //Кл, Na+
    public static double secondIonM = 3.817 * Math.Pow(10, -26); //кг
    
    public Text ui_maxValue;
    public Text ui_minValue;
    public Text ui_midValue;
    
    public Text ui_xPBorderText;
    public Text ui_xNBorderText;
    public Text ui_yPBorderText;
    public Text ui_yNBorderText;
    public Text ui_zPBorderText;
    public Text ui_zNBorderText;
    
    public ParticleSystem ps;
    public GameObject CalcPanel;
    
    public InputField ui_firstIonQ;
    public InputField ui_firstIonM;
    public InputField ui_secondIonQ;
    public InputField ui_secondIonM;

    public List<IonsParameters> ionsParameters = new List<IonsParameters>()
    {
        new IonsParameters(){ q = firstIonQ, m = firstIonM, col = Color.blue}, // Кл, кг
        new IonsParameters(){ q = secondIonQ, m = secondIonM, col = Color.red} // Кл, кг
    };

    public delegate double integratedFunction(double x); //Объявление делегата интегрируемой ф-ций

    private  List<GameObject> linesHolders = new List<GameObject>();

    //=============================================================
    // Инициализация рабочей области
    //=============================================================
    void OnDestroy()
    {
        StateController.ESetState -= OnSetState;
    }

    IEnumerator DestroyPreviousCalculation()
    {
        yield return null;

        foreach (var holder in linesHolders) Destroy(holder);
    }
    
    void Awake()
    {
        StateController.ESetState += OnSetState;
        coilsP = coilsRoot.GetComponentsInChildren<Coil>().ToList();
    }

    private void OnSetState(SceneState state)
    {
        switch (state)
        {
                case SceneState.Setup:
                    StartCoroutine(DestroyPreviousCalculation());
                    break;
                case SceneState.Calculation:
                    StartCoroutine(StartCalculation());
                    break;
        }
    }

    float ConvertFromSIToUnitySize(float value)
    {
        float newValue = value * 10f;
        return newValue;
    }
    
    float ConvertFromUnitySizeToSI(float value)
    {
        float newValue = value / 10f;
        return newValue;
    }
    
    //=============================================================
    // Запуск вычислений по кнопке "включения"
    //=============================================================
    IEnumerator StartCalculation()
    {
        yield return null;

        foreach (var coil in coilsP)
                coil.Reinit();

            //Отрисовываем верные границы области
            drawZone.transform.localScale = new Vector3(xFieldSize * 10f, yFieldSize * 10f, zFieldSize * 10f);
            drawZone.SetActive(true);

            beta = ((doubleAlpha / 2 - Math.PI) / 2);

            // Запуск упрощенного расчета магнитных полей на партиклах
            if (calculatingSimplePowerLines)
            {
                Action CalculateSPL = CalculateSimplifiedPowerLine;
                CalculateSPL();
                // CalculateSPL.BeginInvoke(ar =>
                // {
                //     Debug.LogFormat("Done : {0}", ar.IsCompleted);
                // }, null);
            }
            
            // Запуск расчета магнитных полей
            if (calculatingPowerLines)
            {
                Action CalculatePL = CalculatePowerLine;
                CalculatePL();
                // CalculatePL.BeginInvoke(null, null);
            }
                

            // Запуск расчета траекторий заряженных частиц
            if (calculatingIonTrajectory)
            {
                Action CalculateIT = CalculateIonTrajectory;
                CalculateIT();
                // CalculateIT.BeginInvoke(null, null);
            }   

            SetValuesForDrawZone();
            CalcPanel?.SetActive(false);
    }

    //======================================================
    // Подготовка необходимых заранее вычислений
    //======================================================
    #region
    
    // Нанесение значений на грани рабочей области
    void SetValuesForDrawZone()
    {
        ui_xPBorderText.text = (xFieldSize/2).ToString();
        ui_xNBorderText.text = (-xFieldSize/2).ToString();
        ui_yPBorderText.text = (yFieldSize/2).ToString();
        ui_yNBorderText.text = (-yFieldSize/2).ToString();
        ui_zPBorderText.text = (zFieldSize/2).ToString();
        ui_zNBorderText.text = (-zFieldSize/2).ToString();
        
        // float max3 = Math.Max(xFieldSize, Math.Max(yFieldSize, zFieldSize));

        // ui_xPBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
        // ui_xNBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
        // ui_yPBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
        // ui_yNBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
        // ui_zPBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
        // ui_zNBorderText.transform.localScale = new Vector3(10f / drawZone.transform.localScale.x,
        //     10f / drawZone.transform.localScale.y, 10f / drawZone.transform.localScale.z);
    }

    // Интегрируемая функция K
    double IntegratedFunctionK(double beta)
    {
        return (1 / Math.Sqrt(1 - kSquared(ro, R, z) * Math.Pow(Math.Sin(beta), 2))); 
    }

    // Интегрируемая функция L
    double IntegratedFunctionL(double beta)
    {
        return (Math.Sqrt(1 - kSquared(ro, R, z) * Math.Pow(Math.Sin(beta), 2))); 
    }

    // Переменная для упрощения знаменателя подынтегрального выражения 
    double kSquared(double ro, double R, double z)
    {
        return ((4.0f * ro * R) / (Math.Pow(R + ro, 2) + Math.Pow(z, 2))); 
    }

    void CalculateK()
    {
        // K = TrapezoidalIntegration(IntegratedFunctionK, 0, Math.PI / 2);
        K = SimpsonIntegration(IntegratedFunctionK, 0, Math.PI / 2);
    }

    // Вычисление значения интеграла L
    void CalculateL()
    {
        // L = TrapezoidalIntegration(IntegratedFunctionL, 0, Math.PI / 2);
        L = SimpsonIntegration(IntegratedFunctionL, 0, Math.PI / 2);
    }

    // Вычисление компоненты вектора магнитной индукции Bρ
    double CalculateB_ro()
    {
        return ((mu0 * i / (2 * Math.PI)) * (z / (ro * Math.Sqrt(Math.Pow(R + ro, 2) + Math.Pow(z, 2)))) * (L * (Math.Pow(R, 2) + Math.Pow(ro, 2) + Math.Pow(z, 2)) / (Math.Pow(R - ro, 2) + Math.Pow(z, 2)) - K));
    }
                
    // Вычисление компоненты вектора магнитной индукции Bz
    double CalculateB_z()
    {
        return ((mu0 * i / (2 * Math.PI)) * (1 / (Math.Sqrt(Math.Pow(R + ro, 2) + Math.Pow(z, 2)))) * (L * (Math.Pow(R, 2) - Math.Pow(ro, 2) - Math.Pow(z, 2)) / (Math.Pow(R - ro, 2) + Math.Pow(z, 2)) + K));
    }

    // Интегрирование методом трапеций (аргументы: экземпляр делегата интегрируемой ф-ции, нижний предел интегрирования, верхний предел)
    private double TrapezoidalIntegration(integratedFunction f, double a, double b)
    {
        double I = 0;
        double step = (b - a) / 100;

        double fPrev = f(a);
        double fCurr;
        
        for (double i = a; i < b; i = i + step)
        {
            fCurr = f(i + step);
            I += 0.5 * step * (fPrev + fCurr);
            fPrev = fCurr;
        }

        return I;
    }
    
    // Интегрирование методом Симпсона
    private double SimpsonIntegration(integratedFunction f, double a, double b)
    {
        double I = 0;
        double sum1 = 0d;
        double sum2 = 0d;

        int n = 100;
        double step = (b - a) / n;
        double factor = step/3;
        
        double offset = step;
        int m = 4;
        
        double sum = f(a) + f(b);
        for (int i = 0; i < n - 1; i++)
        {
            sum += m * f(a + offset);
            m = 6 - m;
            offset += step;
        }

        I = factor*sum;
        return I;
    }
    #endregion
    
    //======================================================
    // Проведение основных вычислений потенциала в точке M
    //======================================================
    // Вычисление локальных координат точки поля для катушки (витка катушки)
    // Переписать под МАТРИЧНЫЙ расчёт поворота
    void ToLocalCartesian(ref DQuaternion dQ, Transform pivot)
    {
        Vector3 vec = new Vector3((float)dQ.x, (float)dQ.y, (float)dQ.z);

        // Приводим размеры в Си к Юнити
        Vector3 vecUnity = new Vector3(ConvertFromSIToUnitySize(vec.x), ConvertFromSIToUnitySize(vec.y),
            ConvertFromSIToUnitySize(vec.z));
        
        // Переход
        Vector3 pivotPos = pivot.position;
        vec = vecUnity - pivotPos;
        Vector3 transformVec = pivot.transform.InverseTransformDirection(vec);

        // Приводим размеры Юнити к Си
        Vector3 vecSI = new Vector3(ConvertFromUnitySizeToSI(transformVec.x), ConvertFromUnitySizeToSI(transformVec.y),
            ConvertFromUnitySizeToSI(transformVec.z));
        
        dQ = new DQuaternion(vecSI.x, vecSI.y, vecSI.z, dQ.w, dQ.bx, dQ.by, dQ.bz);
    }

    // Обратное вычисление глобальных координат точки поля для катушки (витка катушки)
    void ToGlobalCartesian(ref DQuaternion dQ, Transform pivot)
    {
        Vector3 vecM = new Vector3((float)dQ.x, (float)dQ.y, (float)dQ.z);
        Vector3 vecB = new Vector3((float)dQ.bx, (float)dQ.by, (float)dQ.bz);

        // Приводим размеры в Си к Юнити
        Vector3 vecMUnity = new Vector3(ConvertFromSIToUnitySize(vecM.x), ConvertFromSIToUnitySize(vecM.y),
            ConvertFromSIToUnitySize(vecM.z));
        Vector3 vecBUnity = new Vector3(ConvertFromSIToUnitySize(vecB.x), ConvertFromSIToUnitySize(vecB.y),
            ConvertFromSIToUnitySize(vecB.z));
        
        // Переход
        Vector3 transformVecM = pivot.transform.TransformDirection(vecMUnity) + pivot.position;
        Vector3 transformVecB = pivot.transform.TransformDirection(vecBUnity);

        // Приводим размеры Юнити к Си
        Vector3 transformvecMSI = new Vector3(ConvertFromUnitySizeToSI(transformVecM.x), ConvertFromUnitySizeToSI(transformVecM.y),
            ConvertFromUnitySizeToSI(transformVecM.z));
        Vector3 transformvecBSI = new Vector3(ConvertFromUnitySizeToSI(transformVecB.x), ConvertFromUnitySizeToSI(transformVecB.y),
            ConvertFromUnitySizeToSI(transformVecB.z));
        
        dQ = new DQuaternion(transformvecMSI.x, transformvecMSI.y, transformvecMSI.z, dQ.w, transformvecBSI.x, transformvecBSI.y, transformvecBSI.z);
    }

    // Переход из локальной ДСК в локальную ЦСК
    void LocalCartesian2LocalCylindrical(DQuaternion M)
    {
        z = M.z;
        ro = Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2));

        if (ro == 0)
            ro = 0.0000001;
    }

    // Производится вычисление в локальной ЦСК компонент Bρ, Bz вектора магнитной индукции в точке M
    void CalculateB()
    {
        CalculateK();
        CalculateL();

        B_ro = CalculateB_ro();
        B_z = CalculateB_z();
    }

    // Переход из локальной ЦСК в локальную ДСК
    DQuaternion LocalCylindrical2LocalCartesian(DQuaternion dQM)
    {
        double den = Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2));
        if (den == 0)
            den = 0.000001;
        B = new DVector3(B_ro * dQM.x / den, B_ro * dQM.y / den, B_z);
        
        // B = new DVector3(B_ro * dQM.x / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))),
        //     B_ro * dQM.y / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))),
        //     B_z);
        DQuaternion newDQ = new DQuaternion(dQM.x, dQM.y, dQM.z, dQM.w, B.x, B.y, B.z);

        return newDQ;
    }

    // Расчёт величины магнитной индукции в точке M
    void CalculateM(ref DQuaternion M)
    {
        M.w = Math.Sqrt(Math.Pow(B.x, 2) + Math.Pow(B.y, 2) + Math.Pow(B.z, 2));
    }

    // // Обсчёт всех значений M в заданном диапазоне
    // void CalculateMagneticField()
    // {
    //     //int index = 0;
    //     maxM = -1000;
    //     minM = 1000;
    //
    //     double localMeshStep = meshStep * 1000f;
    //     double localXFieldSize = xFieldSize * 1000f;
    //
    //     // Запускается обход всех M
    //     for (double xi = -localXFieldSize / 2f; xi <= localXFieldSize / 2f; xi += localMeshStep)
    //     {
    //         //Stopwatch sw = new Stopwatch();
    //         //sw.Start();
    //         //MList.RegisterWorker();
    //         //new Action<double>(DoWorkThreaded).BeginInvoke(xi, null, null);
    //         //Debug.Log("Index : " + index);
    //
    //         //double tempX = xi / 100f;
    //         xi = xi / 100f;
    //         DoWorkThreaded(xi);
    //         xi = xi * 100f;
    //         //DoWorkThreaded(xi);//, ref index);
    //
    //         //sw.Stop();
    //         //Debug.Log("Time elapsed : " + sw.ElapsedMilliseconds);
    //     }
    // }

    // private void DoWorkThreaded(double xi)//, ref int index)
    // {
    //     double localMeshStep = meshStep * 1000f;
    //     double localYFieldSize = yFieldSize * 1000f;
    //     double localZFieldSize = zFieldSize * 1000f;
    //
    //     //int i = index;
    //     for (double yi = -localYFieldSize / 2f; yi <= localYFieldSize / 2f; yi += localMeshStep)
    //     {
    //         for (double zi = -localZFieldSize / 2f; zi <= localZFieldSize / 2f; zi += localMeshStep)
    //         {
    //             MList.RegisterWorker();
    //             //new Action<double, double, double, int>(CalculateInThread).BeginInvoke(xi, yi, zi, i++, null, null);
    //             //new Action<double, double, double>(CalculateInThread).BeginInvoke(xi, yi, zi, null, null);
    //
    //             double tempY = yi / 100f;
    //             double tempZ = zi / 100f;
    //             CalculateInThread(xi, yi / 100f, zi / 100f);
    //             //CalculateInThread(xi, yi, zi);//, i++);
    //
    //             //MList.UnRegisterWorker();
    //         }
    //
    //     }
    //     //index = i;
    // }

    // private void CalculateInThread(double xi, double yi, double zi)//, int index)
    // {
    //     DQuaternion M = new DQuaternion(xi, yi, zi, 0, 0, 0, 0);
    //
    //     LocalCartesian2LocalCylindrical(M);
    //
    //     CalculateB();
    //
    //     M = LocalCylindrical2LocalCartesian(M);
    //
    //     CalculateM(ref M);
    //
    //     if (!double.IsNaN(M.w) || !double.IsInfinity(M.w))
    //     {
    //         // Выкидываются близкие к бесконечности значения
    //         if (M.w < 1000)
    //         {
    //             MList.AddSafely(M);
    //             //MList.AddSafely(M, index);
    //             if ((M.w > maxM) && (!double.IsInfinity(M.w))) maxM = M.w;
    //             if (M.w < minM) minM = M.w;
    //
    //             if (xi > 0)
    //                 maxM = M.w;
    //
    //         }
    //     }
    //     else
    //     {
    //         //Debug.Log("( " + M.x + " , " + M.y + " , " + M.z + " , " + M.w);
    //     }
    //     MList.UnRegisterWorker();
    // }

    // // Настройка визуализации
    // IEnumerator SetupVisualistaionSystem()
    // {
    //     Stopwatch sw = new Stopwatch();
    //     sw.Start();
    //     yield return null;
    //
    //     while (MList.currentWorkersCount > 0)
    //     {
    //         yield return null;
    //     }
    //
    //     int counter = lineArray.Length;
    //
    //     Color lowColor = new Color(0f, 0f, 1f, colorIntensity);
    //     Color maxColor = new Color(1f, 1f, 1f, colorIntensity);
    //
    //     MaterialPropertyBlock props = new MaterialPropertyBlock();
    //
    //     for (int j = 0; j < counter - 1; j++)
    //     {
    //         double localMaxM = maxM - minM;
    //         var tempM = MList.GetPoint(j); //MList.GetPointByIndex(j);
    //
    //         lineArray[j].transform.position = new Vector3((float)tempM.x, (float)tempM.y, (float)tempM.z);
    //         Vector3 orientationPoint = new Vector3((float)tempM.bx * 10000, (float)tempM.by * 10000, (float)tempM.bz * 10000); 
    //         lineArray[j].transform.LookAt(orientationPoint.normalized + lineArray[j].transform.position);
    //
    //         Renderer rnd = lineArray[j].GetComponentInChildren<Renderer>();
    //         props.SetColor("_Color", Color.Lerp(lowColor, maxColor, (float)(2 * tempM.w / localMaxM)));
    //
    //         rnd.SetPropertyBlock(props);
    //         rnd.gameObject.layer = 8;
    //
    //         //Debug.Log("Color: " + prt.color + ", Value: " + (tempM.w / localMaxM));
    //
    //     }
    //
    //     //fieldPoints.SetParticles(particleList.ToArray(), particleList.Count);
    //     sw.Stop();
    //     Debug.LogWarningFormat("Threaded work time seems to : {0}", sw.ElapsedMilliseconds);
    // }

    /// <summary>
    /// Последовательный упрощенный расчет точек линий магнитного поля
    /// </summary>

    private void CalculateSimplifiedPowerLine()
    {
        List<List<LinePoint>> points = new List<List<LinePoint>>();

        int gridStep = 0;

        // int numOfParticles = (int) (((8 * Math.Pow(0.01, 3)) / Math.Pow(localMeshStep, 3)) * 2f);
        int numOfParticles = 5000;
        int pCounter = 0;
        ps.Emit(numOfParticles);

        ParticleSystem.Particle[] pArray = new ParticleSystem.Particle[numOfParticles];
        
        var rings = new List<GameObject>();

        foreach (var coil in coilsP)
            SimplifiedCalculationForParticleSystem(ref gridStep, points, coil, pArray, ref pCounter);

        Debug.Log("Вызов отрисовки");
        
        ps.SetParticles(pArray);
        
        // Material _material = matMag;

        // Debug.Log(points.Count);
        // CreateLinesMesh(points, _material);
    }

    private void SimplifiedCalculationForParticleSystem(ref int gridStep, List<List<LinePoint>> points, Coil growCoil, ParticleSystem.Particle[] pArray, ref int pCounter)
    {
        double localMeshStep = meshStep * 5f;

        gridStep = 0;

        maxM = 0;
        minM = 0;

        // Вычислим значение максимума 
        // Приводим значения к размерам отрисовки Unity
        float val_pos = ConvertFromSIToUnitySize(0.000001f);

        // Делаем отступ на заданный шаг от 0
        var startPos_tmp = new Vector3(val_pos, val_pos, val_pos);
        startPos_tmp = growCoil.transform.TransformDirection(startPos_tmp) + growCoil.transform.position;

        float localXi_tpm = ConvertFromUnitySizeToSI(startPos_tmp.x);
        float localYi_tmp = ConvertFromUnitySizeToSI(startPos_tmp.y);
        float localZi_tmp = ConvertFromUnitySizeToSI(startPos_tmp.z);

        // Подсчёт значения B
        DQuaternion M_i_tmp = new DQuaternion(localXi_tpm, localYi_tmp, localZi_tmp, 0, 0, 0, 0);

        DQuaternion summ_tmp = new DQuaternion();

        foreach (var coil in coilsP)
        {
            foreach (var ring in coil.rings)
            {
                R = coil.R;
                i = coil.i;

                ToLocalCartesian(ref M_i_tmp, ring);
                LocalCartesian2LocalCylindrical(M_i_tmp);
                CalculateB();
                M_i_tmp = LocalCylindrical2LocalCartesian(M_i_tmp);
                CalculateM(ref M_i_tmp);
                ToGlobalCartesian(ref M_i_tmp, ring);

                summ_tmp.x = M_i_tmp.x;
                summ_tmp.y = M_i_tmp.y;
                summ_tmp.z = M_i_tmp.z;

                summ_tmp.w += M_i_tmp.w;
                summ_tmp.bx += M_i_tmp.bx;
                summ_tmp.by += M_i_tmp.by;
                summ_tmp.bz += M_i_tmp.bz;
            }
        }

        // Обновляем значение макс
        maxM = summ_tmp.w;

        // Вычислим значение минимума 
        // Приводим значения к размерам отрисовки Unity
        val_pos = ConvertFromSIToUnitySize((float)growCoil.R * 2);

        // Делаем отступ на заданный шаг от 0
        startPos_tmp = new Vector3(val_pos, val_pos, val_pos);
        startPos_tmp = growCoil.transform.TransformDirection(startPos_tmp) + growCoil.transform.position;

        localXi_tpm = ConvertFromUnitySizeToSI(startPos_tmp.x);
        localYi_tmp = ConvertFromUnitySizeToSI(startPos_tmp.y);
        localZi_tmp = ConvertFromUnitySizeToSI(startPos_tmp.z);

        // Подсчёт значения B
        M_i_tmp = new DQuaternion(localXi_tpm, localYi_tmp, localZi_tmp, 0, 0, 0, 0);

        summ_tmp = new DQuaternion();

        foreach (var coil in coilsP)
        {
            foreach (var ring in coil.rings)
            {
                R = coil.R;
                i = coil.i;

                ToLocalCartesian(ref M_i_tmp, ring);
                LocalCartesian2LocalCylindrical(M_i_tmp);
                CalculateB();
                M_i_tmp = LocalCylindrical2LocalCartesian(M_i_tmp);
                CalculateM(ref M_i_tmp);
                ToGlobalCartesian(ref M_i_tmp, ring);

                summ_tmp.x = M_i_tmp.x;
                summ_tmp.y = M_i_tmp.y;
                summ_tmp.z = M_i_tmp.z;

                summ_tmp.w += M_i_tmp.w;
                summ_tmp.bx += M_i_tmp.bx;
                summ_tmp.by += M_i_tmp.by;
                summ_tmp.bz += M_i_tmp.bz;
            }
        }

        // Обновляем значение макс
        minM = summ_tmp.w;
        
        // Обновляем локальный максимум
        double localMaxM = maxM - minM;
        
        // ПРоизводим вычисления для партиклов
        for (float xi = -(float)growCoil.R * 0.8f; xi <= (float)growCoil.R * 0.8f; xi += (float) localMeshStep) // Окрестность 100 мм, шаг 5 мм
        {
            for (float yi = -(float)growCoil.R * 0.8f; yi <= (float)growCoil.R * 0.8f; yi += (float)localMeshStep)
            {
                for (float zi = -(float)growCoil.turnsNum * (float)growCoil.wireWidth * 0.67f; zi <= (float)growCoil.turnsNum * (float)growCoil.wireWidth * 0.8f; zi += (float) localMeshStep, gridStep++)
                {
                    if (Math.Sqrt(Math.Pow(xi, 2) + Math.Pow(yi, 2)) <= (float)growCoil.R * 0.81f)
                    {
                        // Приводим значения к размерам отрисовки Unity
                        float x_i = ConvertFromSIToUnitySize(xi);
                        float y_i = ConvertFromSIToUnitySize(yi);
                        float z_i = ConvertFromSIToUnitySize(zi);

                        // Делаем отступ на заданный шаг от 0
                        var startPos = new Vector3(x_i, y_i, z_i);
                        startPos = growCoil.transform.TransformDirection(startPos) + growCoil.transform.position;
                        x_i = startPos.x;
                        y_i = startPos.y;
                        z_i = startPos.z;

                        int stepNumber = 0;

                        Color lowColor = new Color(179f / 255f, 222f / 255f, 255f / 255f);
                        Color highColor = new Color(255f / 255f, 73f / 255f, 108f / 255f);
                        lowColor.a = 0.005f;
                        highColor.a = colorIntensity * 0.15f;

                        float xFieldSizeUnity = ConvertFromSIToUnitySize(xFieldSize);
                        float yFieldSizeUnity = ConvertFromSIToUnitySize(yFieldSize);
                        float zFieldSizeUnity = ConvertFromSIToUnitySize(zFieldSize);
                        
                        stepNumber++;

                        float localXi = ConvertFromUnitySizeToSI(x_i);
                        float localYi = ConvertFromUnitySizeToSI(y_i);
                        float localZi = ConvertFromUnitySizeToSI(z_i);

                        // Подсчёт значения B
                        DQuaternion M_i = new DQuaternion(localXi, localYi, localZi, 0, 0, 0, 0);

                        DQuaternion summ = new DQuaternion();

                        foreach (var coil in coilsP)
                        {
                            foreach (var ring in coil.rings)
                            {
                                R = coil.R;
                                i = coil.i;

                                ToLocalCartesian(ref M_i, ring);
                                LocalCartesian2LocalCylindrical(M_i);
                                CalculateB();
                                M_i = LocalCylindrical2LocalCartesian(M_i);
                                CalculateM(ref M_i);
                                ToGlobalCartesian(ref M_i, ring);

                                summ.x = M_i.x;
                                summ.y = M_i.y;
                                summ.z = M_i.z;

                                summ.w += M_i.w;
                                summ.bx += M_i.bx;
                                summ.by += M_i.by;
                                summ.bz += M_i.bz;
                            }
                        }
                        
                        // Визуализация точки
                        var pos = new Vector3(x_i, y_i, z_i);
                        var color = Color.Lerp(lowColor, highColor, (float) (10 * summ.w / localMaxM));
                        // AddLinePoint(points, gridStep, pos, color);
                        // }

                        x_i = ConvertFromSIToUnitySize(localXi);
                        y_i = ConvertFromSIToUnitySize(localYi);
                        z_i = ConvertFromSIToUnitySize(localZi);
                        
                        ParticleSystem.Particle newP = new ParticleSystem.Particle();
                        newP.position = pos;
                        newP.startColor = color;
                        newP.startSize = 0.45f;

                        pArray[pCounter] = newP;
                        pCounter++;
                    }
                }
            }
        }
    }
    
    /// <summary>
    /// Последовательный полный расчет точек линий магнитного поля
    /// </summary>
    private void CalculatePowerLine()
    {
        List<List<LinePoint>> points = new List<List<LinePoint>>();

        int gridStep = 0;

        var rings = new List<GameObject>();

        foreach (var coil in coilsP)
        {
            PointwiseOffset(true, ref gridStep, points, coil);
            PointwiseOffset(false, ref gridStep, points, coil);   
        }

        Material _material = matMag;

        Debug.Log(points.Count);
        CreateLinesMesh(points, _material);
    }
    
    private void PointwiseOffset(bool dir, ref int gridStep, List<List<LinePoint>> points, Coil growCoil)
    {
        double localMeshStep = meshStep;

        maxM = 0;
        minM = 0;

        for (float xi = -0.01f; xi <= 0.01f; xi += (float) localMeshStep) // Окрестность 100 мм, шаг 5 мм
        {
            for (float yi = -0.01f; yi <= 0.01f; yi += (float)localMeshStep, gridStep++)
            {
                // Вырезание окрестностей ноля
                //if ((xi > -4.15f) && ((xi < -3.85f)))
                //{
                //    xi = -xi;
                //    Debug.Log(xi);
                //}
                //if ((yi > -4.15f) && ((yi < -3.85f))) yi = -yi;

                // Приводим значения к размерам отрисовки Unity
                float x_i = ConvertFromSIToUnitySize(xi);
                float y_i = ConvertFromSIToUnitySize(yi);
                float z_i = 0;

                // Делаем отступ на заданный шаг от 0
                var startPos = new Vector3(x_i, y_i, z_i);
                startPos = growCoil.transform.TransformDirection(startPos) + growCoil.transform.position;
                x_i = startPos.x;
                y_i = startPos.y;
                z_i = startPos.z;
                
                int stepNumber = 0;

                Color lowColor = new Color(179f / 255f, 222f / 255f, 255f / 255f);
                Color highColor = new Color(255f / 255f, 73f / 255f, 108f / 255f);
                lowColor.a = colorIntensity;
                highColor.a = colorIntensity;

                float xFieldSizeUnity = ConvertFromSIToUnitySize(xFieldSize);
                float yFieldSizeUnity = ConvertFromSIToUnitySize(yFieldSize);
                float zFieldSizeUnity = ConvertFromSIToUnitySize(zFieldSize);   

                while ((stepNumber < maxSteps) 
                       && (x_i >= (-xFieldSizeUnity / 2f)) && (x_i <= (xFieldSizeUnity / 2f)) 
                       && (y_i >= (-yFieldSizeUnity / 2f)) && (y_i <= (yFieldSizeUnity / 2f)) 
                       && (z_i >= (-zFieldSizeUnity / 2f)) && (z_i <= (zFieldSizeUnity / 2f)))
                {
                    stepNumber++;

                    float localXi = ConvertFromUnitySizeToSI(x_i);
                    float localYi = ConvertFromUnitySizeToSI(y_i);
                    float localZi = ConvertFromUnitySizeToSI(z_i);
                    
                    // Подсчёт значения B
                    DQuaternion M_i = new DQuaternion(localXi, localYi, localZi, 0, 0, 0, 0);

                    DQuaternion summ = new DQuaternion();

                    foreach (var coil in coilsP)
                    {
                        for (int j = 0; j < coil.turnsNum; j++)
                        {
                            foreach (var ring in coil.rings)
                            {
                                R = coil.R - j * coil.wireWidth;
                                i = coil.i;

                                ToLocalCartesian(ref M_i, ring);
                                LocalCartesian2LocalCylindrical(M_i);
                                CalculateB();
                                M_i = LocalCylindrical2LocalCartesian(M_i);
                                CalculateM(ref M_i);
                                ToGlobalCartesian(ref M_i, ring);

                                summ.x = M_i.x;
                                summ.y = M_i.y;
                                summ.z = M_i.z;

                                summ.w += M_i.w;
                                summ.bx += M_i.bx;
                                summ.by += M_i.by;
                                summ.bz += M_i.bz;
                            } 
                        }
                    }

                    // Первый раз обновляем значение макс и мин
                    if (stepNumber == 1)
                    {
                        maxM = summ.w;
                        minM = summ.w;
                    }
                    
                    // Проверка на бесконечные значения
                    if (!double.IsInfinity(summ.w)) Debug.Log("Infinite value on step " + stepNumber + ", point (" + x_i + ", " + y_i + ", " + z_i + ")");
                    
                    if ((summ.w > maxM) && (!double.IsInfinity(summ.w))) maxM = summ.w;
                    if (summ.w < minM) minM = summ.w;

                    double localMaxM = maxM - minM;

                    // Визуализация шага всторону
                    var pos = new Vector3(x_i, y_i, z_i);
                    var color = Color.Lerp(lowColor, highColor, (float) (10 * summ.w / localMaxM));
                    AddLinePoint(points, gridStep, pos, color);
                    
                    // Расчет точек смещения
                    double alpha_i = Math.Sqrt(Math.Pow(summ.bx, 2) + Math.Pow(summ.by, 2) + Math.Pow(summ.bz, 2));
                    
                    // Уменьшение шага смещения от точки к точке
                    float h = 0.01f;
                    localXi = localXi + (h * (float) summ.bx / (float) alpha_i) * (dir ? 1 : -1);
                    localYi = localYi + (h * (float) summ.by / (float) alpha_i) * (dir ? 1 : -1);
                    localZi = localZi + (h * (float) summ.bz / (float) alpha_i) * (dir ? 1 : -1);

                    x_i = ConvertFromSIToUnitySize(localXi);
                    y_i = ConvertFromSIToUnitySize(localYi);
                    z_i = ConvertFromSIToUnitySize(localZi);
                }
            }
        }
        
        // Записываем значения макс и мин значений B в интерфейс
        
        ui_maxValue.text = maxM.ToString("e");
        ui_minValue.text = minM.ToString("e");
        ui_midValue.text = ((maxM + minM)/2).ToString("e");
    }
    
    /// <summary>
    /// Последовательный расчет точек движения заряженной частицы
    /// </summary>
    private void CalculateIonTrajectory()
    {
        DateTime t1 = DateTime.Now;
        
        for (int ii = 0; ii < ionsParameters.Count; ii++)
        {
            List<List<LinePoint>> ionPoints = new List<List<LinePoint>>();

            double x_prev = p0.x; // м
            double y_prev = p0.y; // м
            double z_prev = p0.z; // м
            double x = p1.x; // м
            double y = p1.y; // м
            double z = p1.z; // м

            // h = Math.Pow(10, -6);
            h = 2 * Math.Pow(10, -7);

            for (int j = 0; j < ionSteps; j++)
            {
                // float localXi = ConvertFromUnitySizeToSI(x);
                // float localYi = ConvertFromUnitySizeToSI(y);
                // float localZi = ConvertFromUnitySizeToSI(z);
                    
                // Подсчёт значения B
                DQuaternion M_i = new DQuaternion(x, y, z, 0, 0, 0, 0);

                DQuaternion summ = new DQuaternion();

                foreach (var coil in coilsP)
                {
                    for (int jj = 0; jj < coil.turnsNum; jj++)
                    {
                        foreach (var ring in coil.rings)
                        {
                            R = coil.R - jj * coil.wireWidth;
                            i = coil.i;

                            ToLocalCartesian(ref M_i, ring);
                            LocalCartesian2LocalCylindrical(M_i);
                            CalculateB();
                            M_i = LocalCylindrical2LocalCartesian(M_i);
                            CalculateM(ref M_i);
                            ToGlobalCartesian(ref M_i, ring);

                            summ.x = M_i.x;
                            summ.y = M_i.y;
                            summ.z = M_i.z;

                            summ.w += M_i.w;
                            summ.bx += M_i.bx;
                            summ.by += M_i.by;
                            summ.bz += M_i.bz;
                        } 
                    }
                }

                // Расчет электрической составляющей
                Vector3 EVector = new Vector3(-1f, 0.5f, 0.5f); // Задаём единичным вектором
                double E = 0.00001d; // В/м
                // double E = 1d / 1000; // В/м
                
                double alphaTemp = ionsParameters[ii].q * h / (2 * ionsParameters[ii].m);
                
                double x_next = 0;
                double y_next = 0;
                double z_next = 0;

                // Случай с магнитным полем, без электрического
                if (mField)
                {
                    double c1 = 2 * x - x_prev - alphaTemp * summ.bz * y_prev + alphaTemp * summ.by * z_prev;
                    double c2 = alphaTemp * summ.bx * x_prev + 2 * y - y_prev - alphaTemp * summ.bx * z_prev;
                    double c3 = - alphaTemp * summ.by * x_prev + alphaTemp * summ.bx * y_prev + 2 * z - z_prev;
                    
                    double fDelta = 1 + Math.Pow(alphaTemp, 2) * (Math.Pow(summ.bx, 2) + Math.Pow(summ.by, 2) + Math.Pow(summ.bz, 2));
                    
                    double deltaU1 = c1 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bx, 2)) + c2 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (alphaTemp * summ.bx * summ.bz - alphaTemp * Math.Pow(summ.by, 2));
                    double deltaU2 = - c1 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.bx * summ.by) + c2 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.by, 2)) + c3 * (alphaTemp * summ.bx + Math.Pow(alphaTemp, 2) * summ.by * summ.bz);
                    double deltaU3 = c1 * (Math.Pow(alphaTemp, 2) * summ.bx * summ.bz + alphaTemp * summ.by) - c2 * (alphaTemp * summ.bx - Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bz, 2));
                    
                    x_next = deltaU1 / fDelta;
                    y_next = deltaU2 / fDelta;
                    z_next = deltaU3 / fDelta;
                }

                // Случай с магнитным и постоянным электрическим полями
                if (mPlusEFields)
                {
                    double c1 = 2 * x - x_prev - alphaTemp * summ.bz * y_prev + alphaTemp * summ.by * z_prev;
                    double c2 = alphaTemp * summ.bx * x_prev + 2 * y - y_prev - alphaTemp * summ.bx * z_prev;
                    double c3 = - alphaTemp * summ.by * x_prev + alphaTemp * summ.bx * y_prev + 2 * z - z_prev;
                    
                    double fDelta = 1 + Math.Pow(alphaTemp, 2) * (Math.Pow(summ.bx, 2) + Math.Pow(summ.by, 2) + Math.Pow(summ.bz, 2));
                    
                    // Добавить + E_x(xi, yi, zi), + E_y(xi, yi, zi), + E_z(xi, yi, zi)
                    double deltaU1 = c1 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bx, 2)) + c2 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (alphaTemp * summ.bx * summ.bz - alphaTemp * Math.Pow(summ.by, 2)) + EVector.x * E;
                    double deltaU2 = - c1 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.bx * summ.by) + c2 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.by, 2)) + c3 * (alphaTemp * summ.bx + Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + EVector.y * E;
                    double deltaU3 = c1 * (Math.Pow(alphaTemp, 2) * summ.bx * summ.bz + alphaTemp * summ.by) - c2 * (alphaTemp * summ.bx - Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bz, 2)) + EVector.z * E;
                    
                    x_next = deltaU1 / fDelta;
                    y_next = deltaU2 / fDelta;
                    z_next = deltaU3 / fDelta;
                }

                // // Случай с магнитным и переменным электрическим полями
                // if (mPlusDEFields)
                // {
                //     double fDelta = 1 + Math.Pow(K_i_elmag, 2);
                //     double delta_U1 = 2 * x - x_prev - K_i_elmag * y_prev + L_i_elmag * sinGamma * cosBeta + K_i_elmag * (2 * y - y_prev + K_i_elmag * x_prev + L_i_elmag * sinGamma * sinBeta);
                //     double delta_U2 = 2 * y - y_prev + K_i_elmag * x_prev + L_i_elmag * sinGamma * sinBeta - K_i_elmag * (2 * x - x_prev - K_i_elmag * y_prev + L_i_elmag * sinGamma * cosBeta);
                //     double delta_U3 = 2 * y - y_prev + K_i_elmag * x_prev + L_i_elmag * sinGamma * sinBeta - K_i_elmag * (2 * x - x_prev - K_i_elmag * y_prev + L_i_elmag * sinGamma * cosBeta);
                //     x_next = delta_U1 / fDelta;
                //     y_next = delta_U2 / fDelta;
                //     z_next = delta_U3 / fDelta;
                // }
                
                // Вчисление расстояния до самой удаленной точки для масштабирования
                // double maxOfThree = Math.Max(x, Math.Max(y, z));
                // if (maxOfThree > maxIonRange)
                //     maxIonRange = (float)maxOfThree;

                // Приводим размеры в Си к Юнити
                float xUnitySize = ConvertFromSIToUnitySize((float)x);
                float yUnitySize = ConvertFromSIToUnitySize((float)y);
                float zUnitySize = ConvertFromSIToUnitySize((float)z);

                //Рисование
                Vector3 pos = new Vector3(xUnitySize, yUnitySize, zUnitySize);
                AddLinePoint(ionPoints, 1, pos, ionsParameters[ii].col);

                Debug.Log("x diff = " + (x_next - x) + ", y diff = " + (y_next - y) + ", z diff = " + (z_next - z));

                // Переприсвоение перед новым циклом обсчёта
                x_prev = x; 
                y_prev = y; 
                z_prev = z;
                x = x_next; 
                y = y_next; 
                z = z_next;
            }

            Material _material;

            switch (ii)
            {
                case 1:
                    _material = matLine3;
                    break;
                case 2:
                    _material = matLine2;
                    break;
                default:
                    _material = matLine1;
                    break;
            }

            Debug.Log(ionPoints.Count);
            CreateLinesMesh(ionPoints, _material);

            //CreateCompressedLinesMesh(ionPoints, ref meshHolderIon);
        }

        DateTime t2 = DateTime.Now;
        string s = (t2 - t1).ToString();
        Debug.Log(s);
    }

    private void AddLinePoint(List<List<LinePoint>> list, int step, Vector3 pos, Color color)
    {
        if (list.Count < step + 1)
            while(list.Count < step+1)
            {list.Add(new List<LinePoint>());}
        list[step].Add(new LinePoint(){pos = pos, color = color});
    }

    public struct LinePoint
    {
        public Vector3 pos;
        public Color color;
    }
    
    private void CreateLinesMesh(List<List<LinePoint>> points, Material mtr)
    {
        var linesHolder = new GameObject("LinesHolder");
        linesHolders.Add(linesHolder);

        for (int i = points.Count - 1; i >= 0; i--)
        {
            if (points[i].Count < 2)
            {
                points.RemoveAt(i);
                continue;
            }
        }

        for (var x = 0; x < points.Count; x++)
        {
            List<Vector3> positions = new List<Vector3>();
            List<Vector3> normals = new List<Vector3>();
            List<Vector2> uvs = new List<Vector2>();
            List<Color> colors = new List<Color>();
            List<int> indexes = new List<int>();

            for (int y = 0; y < points[x].Count; y++)
            {
                positions.Add(points[x][y].pos);
                positions.Add(points[x][y].pos);

                colors.Add(points[x][y].color);
                colors.Add(points[x][y].color);

                uvs.Add(Vector2.left);
                uvs.Add(-Vector2.left);
            }

            for (int ind = 0; ind < points[x].Count - 1; ind++)
            {
                indexes.Add(ind * 2);
                indexes.Add(ind * 2 + 1);
                indexes.Add(ind * 2 + 3);
                indexes.Add(ind * 2 + 2);
            }
            
            for (int ind = 0; ind < points[x].Count - 1; ind++)
            {
                var nrm = points[x][ind + 1].pos - points[x][ind].pos;
                normals.Add(nrm.normalized);
                normals.Add(nrm.normalized);
            }

            var n = (points[x][points[x].Count - 1].pos - points[x][points[x].Count - 2].pos).normalized;
            normals.Add(n);
            normals.Add(n);

            var mesh = new Mesh();
            mesh.bounds = new Bounds(Vector3.zero, Vector3.one * 100500);
            mesh.indexFormat = IndexFormat.UInt32;
            mesh.vertices = positions.ToArray();
            mesh.normals = normals.ToArray();
            mesh.uv = uvs.ToArray();
            mesh.colors = colors.ToArray();
            mesh.SetIndices(indexes.ToArray(), MeshTopology.Quads, 0, false);

            var meshHolderObj = new GameObject("MeshHolder");
            meshHolderObj.transform.SetParent(linesHolder.transform);
            var filter = meshHolderObj.AddComponent<MeshFilter>();
            filter.sharedMesh = mesh;
            var renderer = meshHolderObj.AddComponent<MeshRenderer>();
            mtr.color = Color.gray;
            renderer.sharedMaterial = mtr;
            renderer.receiveShadows = false;
            renderer.shadowCastingMode = ShadowCastingMode.Off;
        }
    }

    ///////////
    /// Старый кусок кода для отрисовки элементов кривых
    ///////////
    // private void SpawnGoSegment(Vector3 pos, DQuaternion summ, Color color)
    // {
    //     GameObject line_i_right = Instantiate(lineSource);
    //
    //     line_i_right.transform.parent = fieldSource.transform;
    //     line_i_right.transform.position = pos;
    //     
    //     Vector3 orientationPointRight =
    //         new Vector3((float) summ.bx * 10000, (float) summ.@by * 10000, (float) summ.bz * 10000);
    //     
    //     line_i_right.transform.LookAt(orientationPointRight.normalized + line_i_right.transform.position);
    //     
    //     MaterialPropertyBlock props2 = new MaterialPropertyBlock();
    //     Renderer rnd2 = line_i_right.GetComponentInChildren<Renderer>();
    //
    //     props2.SetColor("_Color", color);
    //     rnd2.SetPropertyBlock(props2);
    //     rnd2.gameObject.layer = 8;
    // }
    
    ///////////
    /// Старый кусок кода для генерации сжатого меша
    ///////////
    // private void CreateCompressedLinesMesh(List<List<LinePoint>> points, ref GameObject meshHolderObj)
    // {
    //     for (int i = points.Count - 1; i >= 0; i--)
    //     {
    //         if (points[i].Count < 2)
    //         {
    //             points.RemoveAt(i);
    //             continue;
    //         }
    //     }
    //
    //     List<Vector3> positions = new List<Vector3>();
    //     List<int> indexes = new List<int>();
    //     List<Color> colors = new List<Color>();
    //
    //     var indexCounter = 0;
    //
    //     float scaleCoeff = 1f;
    //
    //     for (int x = 0; x < points.Count; x++)
    //     {
    //         for (int y = 1; y < points[x].Count; y++)
    //         {
    //
    //             Vector3 tempVecPos1 = new Vector3(points[x][y - 1].pos.x * scaleCoeff, points[x][y - 1].pos.y * scaleCoeff, points[x][y - 1].pos.z * scaleCoeff);
    //             Vector3 tempVecPos2 = new Vector3(points[x][y].pos.x * scaleCoeff, points[x][y].pos.y * scaleCoeff, points[x][y].pos.z * scaleCoeff);
    //
    //             positions.Add(tempVecPos1);
    //             positions.Add(tempVecPos2);
    //             colors.Add(points[x][y - 1].color);
    //             colors.Add(points[x][y].color);
    //             indexes.Add(indexCounter++);
    //             indexes.Add(indexCounter++);
    //         }
    //     }
    //
    //     var mesh = new Mesh();
    //     mesh.bounds = new Bounds(Vector3.zero, Vector3.one * 100500);
    //     mesh.indexFormat = IndexFormat.UInt32;
    //     mesh.vertices = positions.ToArray();
    //     mesh.colors = colors.ToArray();
    //     mesh.SetIndices(indexes.ToArray(), MeshTopology.Lines, 0, false);
    //
    //     meshHolderObj = new GameObject("MeshHolder");
    //     var filter = meshHolderObj.AddComponent<MeshFilter>();
    //     
    //     filter.sharedMesh = mesh;
    //     
    //     var renderer = meshHolderObj.AddComponent<MeshRenderer>();
    //     var shader = Shader.Find("Particles/Standard Unlit");
    //     var mat = new Material(shader);
    //     
    //     mat.color = Color.gray;
    //     
    //     renderer.sharedMaterial = mat;
    //     renderer.receiveShadows = false;
    //     renderer.shadowCastingMode = ShadowCastingMode.Off;
    // }
}
