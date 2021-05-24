using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.Linq;
using System.Threading;
using UnityEngine;
using UnityEngine.Rendering;
using UnityEngine.UI;
using Debug = UnityEngine.Debug;
using Random = UnityEngine.Random;

public class IonsParams {
    public double q; 
    public double m;
    public Color col;
}

public class PointMagneticInduction : MonoBehaviour
{
    private static CultureInfo Culture = CultureInfo.InvariantCulture;
    
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

    private double z; //
    private double R; // радиус контура S
    private double i; // величина силы тока, протекающего по контуру
    public double O = 0.0f; // центр контура S
    private double ro; // расстояние от точки M до оси Oz
    public double dl = 0.0f; // элементарный участок контура
    public double doubleAlpha; // угол между двумя элементарными участками контура, расположенными симметрично
    private double beta;
    private double r; // расстояние от точки M до элементарных участков контура dl
    public double mu0; // магнитная постоянная

    private double K;
    private double L;
    private double B_ro;
    private double B_z;

    public double maxM;
    public double minM;

    public double meshStep; // шаг обсчёта
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

    public string ionName;
    public int ionSteps;
    public InputField ui_ionSteps;
    
    public Vector3 p0;
    public InputField p0_X, p0_Y, p0_Z;
    public Vector3 p1;
    public InputField p1_X, p1_Y, p1_Z;
    
    public float gammaAngle;
    public float betaAngle;
    float maxIonRange = 0;
    
    public static double firstIonQ = -1.602 * Math.Pow(10, -19);
    public static double firstIonM = 9.1 * Math.Pow(10, -31);
    public static double secondIonQ = -1.602 * Math.Pow(10, -19);
    public static double secondIonM = 9.1 * Math.Pow(10, -31);
    
    public InputField ui_firstIonQ;
    public InputField ui_firstIonM;
    public InputField ui_secondIonQ;
    public InputField ui_secondIonM;

    public List<IonsParams> ionsParams = new List<IonsParams>()
    {
        new IonsParams(){ q = firstIonQ, m = firstIonM, col = Color.blue}, // Кл, кг
        //new IonsParams(){ q = -2 * 1.602 * Math.Pow(10, -19), m = 1.66 * Math.Pow(10, -27), col = Color.red},
        new IonsParams(){ q = secondIonQ, m = secondIonM, col = Color.red}//,
        //new IonsParams(){ q = 3.817 * Math.Pow(10, -26), m = 3.82 * Math.Pow(10, -23), col = Color.magenta}//,
        //new IonsParams(){ q = -1.602 * Math.Pow(10, -19), m = 9.1 * Math.Pow(10, -31), col = Color.blue}
    };

    public delegate double integratedFunction(double x); //Объявление делегата интегрируемой ф-ций

    private  List<GameObject> linesHolders = new List<GameObject>();

    void Start()
    {
        GetUIValues();
    }

    public void GetUIValues()
    {
        StartCoroutine("GetUiValuesLater");
    }

    IEnumerator GetUiValuesLater()
    {
        yield return new WaitForEndOfFrame();
        if (!(powerLinesToggle is null)) calculatingPowerLines = powerLinesToggle.isOn;
        if (!(ionTrajectoryToggle is null)) calculatingIonTrajectory = ionTrajectoryToggle.isOn;
        if (!(mFieldToggle is null)) mField = mFieldToggle.isOn;
        if (!(mPlusEFieldsToggle is null)) mPlusEFields = mPlusEFieldsToggle.isOn;
        if (!(mPlusDEFieldsToggle is null)) mPlusDEFields = mPlusDEFieldsToggle.isOn;

        if (!(p0_X is null)) p0.x = Convert.ToSingle(p0_X.text, Culture);
        if (!(p0_Y is null)) p0.y = Convert.ToSingle(p0_Y.text, Culture);
        if (!(p0_Z is null)) p0.z = Convert.ToSingle(p0_Z.text, Culture);
        
        if (!(p1_X is null)) p1.x = Convert.ToSingle(p1_X.text, Culture);
        if (!(p1_Y is null)) p1.y = Convert.ToSingle(p1_Y.text, Culture);
        if (!(p1_Z is null)) p1.z = Convert.ToSingle(p1_Z.text, Culture);
        
        if (!(ui_xFieldSize is null))xFieldSize = Convert.ToSingle(ui_xFieldSize.text, Culture);
        if (!(ui_yFieldSize is null))yFieldSize = Convert.ToSingle(ui_yFieldSize.text, Culture);
        if (!(ui_zFieldSize is null))zFieldSize = Convert.ToSingle(ui_zFieldSize.text, Culture);

        Double mesh_StepSizeL;
        if (!(ui_meshStep is null) && Double.TryParse(ui_meshStep.text, out mesh_StepSizeL)) meshStep = mesh_StepSizeL; // А зачем тут double переменная? Может и float хватит?
        if (!(ui_maxSteps is null))maxSteps = Convert.ToInt32(ui_maxSteps.text, Culture);
        if (!(ui_maxSteps is null))ionSteps = Convert.ToInt32(ui_ionSteps.text, Culture);
        
        if (!(ui_firstIonQ is null) && Double.TryParse(ui_firstIonQ.text, out mesh_StepSizeL)) firstIonQ = mesh_StepSizeL; // А зачем тут double переменная? Может и float хватит?
        if (!(ui_firstIonM is null) && Double.TryParse(ui_firstIonM.text, out mesh_StepSizeL)) firstIonM = mesh_StepSizeL;
        if (!(ui_secondIonQ is null) && Double.TryParse(ui_secondIonQ.text, out mesh_StepSizeL)) secondIonQ = mesh_StepSizeL;
        if (!(ui_secondIonM is null) && Double.TryParse(ui_secondIonM.text, out mesh_StepSizeL)) secondIonM = mesh_StepSizeL;
    }
    
    // Use this for initialization
    void OnDestroy()
    {
        StateController.ESetState -= OnSetState;
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
                    StartCoroutine(DestroyPrevCalculation());
                    break;
                case SceneState.Calculation:
                    StartCoroutine(StartCalculation());
                    break;
        }
    }

    IEnumerator DestroyPrevCalculation()
    {
        yield return null;

        foreach (var holder in linesHolders) Destroy(holder);
    }

    //=============================================================
    // Запуск вычислений по кнопке "включения"
    //=============================================================
    IEnumerator StartCalculation()
    {
        yield return null;

        foreach (var coil in coilsP)
            coil.Reinit();

        Stopwatch sw = new Stopwatch();
        sw.Start();

        //Отрисовываем верные границы области
        drawZone.transform.localScale = new Vector3(xFieldSize * 10f, yFieldSize * 10f, zFieldSize * 10f);
        drawZone.SetActive(true);

        mu0 = (4 * Math.PI * Math.Pow(10, -7));

        beta = ((doubleAlpha / 2 - Math.PI) / 2);

        // БЫЛО ВКЛЮЧЕНО (17.09)
        //int workerThreads, completionPortThreads;
        //ThreadPool.GetMaxThreads(out workerThreads, out completionPortThreads);
        //Debug.LogWarningFormat("Доступно потоков : {0}", workerThreads);
        ////new Action(CalculateMagneticField).BeginInvoke(null, null);
        //CalculateMagneticField();

        // Отрисовка
        //InstantiateVisualisationObjects();
        //StartCoroutine(SetupVisualistaionSystem());

        Vector3 vector = new Vector3(1, 2, 3);

        // Запуск расчета магнитных полей
        if (calculatingPowerLines)
            new Action(CalculatePowerLine)();//.BeginInvoke(null, null);

        // Запуск расчета траекторий заряженных частиц
        if (calculatingIonTrajectory)
            new Action(CalculateIonTrajectory)();//.BeginInvoke(Callback, null);

        //StartCoroutine(SetupParticleSystem());
        sw.Stop();
        //Debug.LogWarningFormat("Time to calculate : {0}", sw.ElapsedMilliseconds);

        //StateController.SetState(SceneState.Spectate);
    }

    private void Callback(IAsyncResult ar)
    {
        Debug.Log("----------------DONE");
        StateController.SetState(SceneState.Spectate);
    }

    //=============================================================
    // Подготовка рабочей области и заполнение массива точек
    //=============================================================
    //void instantiatePointField()
    //{
    //    for (float x = -xFieldSize / 2f; x <= xFieldSize / 2f; x += (float)meshStep)
    //        for (float y = -yFieldSize / 2f; y <= yFieldSize / 2f; y += (float)meshStep)
    //            for (float z = -zFieldSize / 2f; z <= zFieldSize / 2f; z += (float)meshStep)
    //            {
    //                DQuaternion M = new DQuaternion(x, y, z, 0, 0, 0, 0);
    //                MList.AddSafely(M);
    //            }
    //}

    //======================================================
    // Подготовка необходимых заранее вычислений
    //======================================================
    #region

    // Интегрируемая функция K
    double integratedFunctionK(double beta)
    {
        return (1 / Math.Sqrt(1 - kSquared(ro, R, z) * Math.Pow(Math.Sin(beta), 2))); 
    }

    // Интегрируемая функция L
    double integratedFunctionL(double beta)
    {
        return (Math.Sqrt(1 - kSquared(ro, R, z) * Math.Pow(Math.Sin(beta), 2))); 
    }

    // Переменная для упрощения знаменателя подынтегрального выражения 
    double kSquared(double ro, double R, double z)
    {
        return ((4.0 * ro * R) / (Math.Pow(R + ro, 2) + Math.Pow(z, 2))); 
    }

    void CalculateK()
    {
        K = integration(integratedFunctionK, 0, Math.PI / 2);
    }

    // Вычисление значения интеграла L
    void CalculateL()
    {
        L = integration(integratedFunctionL, 0, Math.PI / 2);
    }

    //// Вычисление векторного потенциала A в точке M
    //double CalculateA()
    //{
    //    return ((mu0 * i / (2 * Math.PI)) * (Math.Sqrt(R / ro)) * ((2 - kSquared(ro, R, z) * K / Math.Sqrt(kSquared(ro, R, z)) - (2 * L) / Math.Sqrt(kSquared(ro, R, z)))));
    //}

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
    public double integration(integratedFunction f, double a, double b)
    {
        double I = 0;
        double step = (b - a) / 100;

        for (double i = a; i < b; i = i + step)
        {
            I = I + ((f(i) + f(i + step)) / 2) * ((i + step) - i);
        }

        return I;
    }
    #endregion

    //======================================================
    // Инстансинг объектов отрисовки
    //======================================================
    //void InstantiateVisualisationObjects()
    //{
    //    int lineCounter = (int)(xFieldSize/meshStep + 1) * (int)(yFieldSize / meshStep + 1) * (int)(zFieldSize / meshStep + 1);
    //    lineArray = new GameObject[lineCounter];

    //    for (int i = 0; i < lineCounter - 1; i++)
    //    {
    //        lineArray[i] = Instantiate(lineSource);
    //        lineArray[i].transform.parent = fieldSource.transform;
    //    }
    //}

    //======================================================
    // Проведение основных вычислений потенциала в точке M
    //======================================================
    // Вычисление локальных координат точки поля для катушки (витка катушки)
    // Переписать под ВЕКТОРНЫЙ расчёт поворота
    void ToLocalCartesian(ref DQuaternion dQ, Transform pivot)
    {
        Vector3 vec = new Vector3((float)dQ.x, (float)dQ.y, (float)dQ.z);
        vec = vec - pivot.position;
        Vector3 transformVec = pivot.transform.InverseTransformDirection(vec);

        dQ = new DQuaternion(transformVec.x, transformVec.y, transformVec.z, dQ.w, dQ.bx, dQ.by, dQ.bz);
    }

    // Обратное вычисление глобальных координат точки поля для катушки (витка катушки)
    void ToGlobalCartesian(ref DQuaternion dQ, Transform pivot)
    {
        Vector3 vecM = new Vector3((float)dQ.x, (float)dQ.y, (float)dQ.z);
        Vector3 vecB = new Vector3((float)dQ.bx, (float)dQ.by, (float)dQ.bz);

        Vector3 transformVecM = pivot.transform.TransformDirection(vecM) + pivot.position;
        Vector3 transformVecB = pivot.transform.TransformDirection(vecB);

        dQ = new DQuaternion(transformVecM.x, transformVecM.y, transformVecM.z, dQ.w, transformVecB.x, transformVecB.y, transformVecB.z);
    }

    // Переход из локальной ДСК в локальную ЦСК
    void LocalCartesian2LocalCylindrical(DQuaternion M)
    {
        z = M.z;
        ro = Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2));
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
    //void LocalCylindrical2LocalCartesian(ref DQuaternion M)
    //{
    //    B = new DVector3(B_ro * M.x / (Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2))), B_ro * M.y / (Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2))), B_z);

    //    M.bx = B.x;
    //    M.by = B.y;
    //    M.bz = B.z;
    //}

    //Перегрузка перехода из локальной ЦСК в локальную ДСК
    DQuaternion LocalCylindrical2LocalCartesian(DQuaternion dQM)
    {
        B = new DVector3(B_ro * dQM.x / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))),
            B_ro * dQM.y / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))),
            B_z);
        DQuaternion newDQ = new DQuaternion(dQM.x, dQM.y, dQM.z, dQM.w, B.x, B.y, B.z);

        return newDQ;
    }

    // Расчёт величины магнитной индукции в точке M
    void CalculateM(ref DQuaternion M)
    {
        M.w = Math.Sqrt(Math.Pow(B.x, 2) + Math.Pow(B.y, 2) + Math.Pow(B.z, 2));
    }

    // Обсчёт всех значений M в заданном диапазоне
    void CalculateMagneticField()
    {
        //int index = 0;
        maxM = -1000;
        minM = 1000;

        double localMeshStep = meshStep * 1000f;
        double localXFieldSize = xFieldSize * 1000f;

        // Запускается обход всех M
        for (double xi = -localXFieldSize / 2f; xi <= localXFieldSize / 2f; xi += localMeshStep)
        {
            //Stopwatch sw = new Stopwatch();
            //sw.Start();
            //MList.RegisterWorker();
            //new Action<double>(DoWorkThreaded).BeginInvoke(xi, null, null);
            //Debug.Log("Index : " + index);

            //double tempX = xi / 100f;
            xi = xi / 100f;
            DoWorkThreaded(xi);
            xi = xi * 100f;
            //DoWorkThreaded(xi);//, ref index);

            //sw.Stop();
            //Debug.Log("Time elapsed : " + sw.ElapsedMilliseconds);
        }
    }

    private void DoWorkThreaded(double xi)//, ref int index)
    {
        double localMeshStep = meshStep * 1000f;
        double localYFieldSize = yFieldSize * 1000f;
        double localZFieldSize = zFieldSize * 1000f;

        //int i = index;
        for (double yi = -localYFieldSize / 2f; yi <= localYFieldSize / 2f; yi += localMeshStep)
        {
            for (double zi = -localZFieldSize / 2f; zi <= localZFieldSize / 2f; zi += localMeshStep)
            {
                MList.RegisterWorker();
                //new Action<double, double, double, int>(CalculateInThread).BeginInvoke(xi, yi, zi, i++, null, null);
                //new Action<double, double, double>(CalculateInThread).BeginInvoke(xi, yi, zi, null, null);

                double tempY = yi / 100f;
                double tempZ = zi / 100f;
                CalculateInThread(xi, yi / 100f, zi / 100f);
                //CalculateInThread(xi, yi, zi);//, i++);

                //MList.UnRegisterWorker();
            }

        }
        //index = i;
    }

    private void CalculateInThread(double xi, double yi, double zi)//, int index)
    {
        DQuaternion M = new DQuaternion(xi, yi, zi, 0, 0, 0, 0);

        LocalCartesian2LocalCylindrical(M);

        CalculateB();

        M = LocalCylindrical2LocalCartesian(M);

        CalculateM(ref M);

        if (!double.IsNaN(M.w) || !double.IsInfinity(M.w))
        {
            // Выкидываются близкие к бесконечности значения
            if (M.w < 1000)
            {
                MList.AddSafely(M);
                //MList.AddSafely(M, index);
                if ((M.w > maxM) && (!double.IsInfinity(M.w))) maxM = M.w;
                if (M.w < minM) minM = M.w;

                if (xi > 0)
                    maxM = M.w;

            }
        }
        else
        {
            //Debug.Log("( " + M.x + " , " + M.y + " , " + M.z + " , " + M.w);
        }
        MList.UnRegisterWorker();
    }

    // Настройка визуализации
    IEnumerator SetupVisualistaionSystem()
    {
        Stopwatch sw = new Stopwatch();
        sw.Start();
        yield return null;

        while (MList.currentWorkersCount > 0)
        {
            yield return null;
        }

        //yield return new WaitForSeconds(1f);

        //while (MList.currentWorkersCount > 0)
        //{
        //    yield return null;
        //}

        int counter = lineArray.Length;

        Color lowColor = new Color(0f, 0f, 1f, colorIntensity);
        Color maxColor = new Color(1f, 1f, 1f, colorIntensity);

        MaterialPropertyBlock props = new MaterialPropertyBlock();

        for (int j = 0; j < counter - 1; j++)
        {
            double localMaxM = maxM - minM;
            var tempM = MList.GetPoint(j); //MList.GetPointByIndex(j);

            lineArray[j].transform.position = new Vector3((float)tempM.x, (float)tempM.y, (float)tempM.z);
            Vector3 orientationPoint = new Vector3((float)tempM.bx * 10000, (float)tempM.by * 10000, (float)tempM.bz * 10000); 
            lineArray[j].transform.LookAt(orientationPoint.normalized + lineArray[j].transform.position);

            Renderer rnd = lineArray[j].GetComponentInChildren<Renderer>();
            props.SetColor("_Color", Color.Lerp(lowColor, maxColor, (float)(2 * tempM.w / localMaxM)));

            rnd.SetPropertyBlock(props);
            rnd.gameObject.layer = 8;

            //Debug.Log("Color: " + prt.color + ", Value: " + (tempM.w / localMaxM));

        }

        //fieldPoints.SetParticles(particleList.ToArray(), particleList.Count);
        sw.Stop();
        Debug.LogWarningFormat("Threaded work time seems to : {0}", sw.ElapsedMilliseconds);
    
    
    }

    // Расчет точек движения заряженной частицы
    private void CalculateIonTrajectory()
    {
        for (int ii = 0; ii < ionsParams.Count; ii++)
        {
            List<List<LinePoint>> ionPoints = new List<List<LinePoint>>();

            double x_prev = p0.x / 1000f; double y_prev = p0.y / 1000f; double z_prev = p0.z / 1000f; //в метрах
            double x = p1.x / 1000f; double y = p1.y / 1000f; double z = p1.z / 1000f; //в метрах

            double h = Math.Pow(10, -6);
            //double h = Math.Pow(10, -4);

            for (int j = 0; j < ionSteps; j++)
            {
                // Подсчёт значения B
                DQuaternion M_i = new DQuaternion(x, y, z, 0, 0, 0, 0);

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

                        summ.w += M_i.w; // * coil.turnsNum; // умножение на количество намоток
                        summ.bx += M_i.bx;
                        summ.by += M_i.by;
                        summ.bz += M_i.bz;
                    }
                }

                // double K_i_mag = (ionsParams[ii].q * h / (2 * ionsParams[ii].m)) * summ.w;
                // double K_i_elmag = (ionsParams[ii].q * h / (2 * ionsParams[ii].m)) * summ.w * Math.Cos(omega * j * h);
                // //double K_p = (C_ / 2 * h) * summ.w * Math.Cos(omega3); // Не хватает домножения на t

                // Расчет электрической составляющей
                //double E = 1; // В/м
                double E = 1d / 1000; // В/м
                //double L_i = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E * Math.Cos(omega * j * h);
                // double L_i_mag = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E;
                // double L_i_elmag = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E * Math.Cos(omega * j * h);

                double alphaTemp = ionsParams[ii].q * h / (2 * ionsParams[ii].m);
                
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
                    
                    double delta_U1 = c1 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bx, 2)) + c2 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (alphaTemp * summ.bx * summ.bz - alphaTemp * Math.Pow(summ.by, 2));
                    double delta_U2 = - c1 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.bx * summ.by) + c2 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.by, 2)) + c3 * (alphaTemp * summ.bx + Math.Pow(alphaTemp, 2) * summ.by * summ.bz);
                    double delta_U3 = c1 * (Math.Pow(alphaTemp, 2) * summ.bx * summ.bz + alphaTemp * summ.by) - c2 * (alphaTemp * summ.bx - Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bz, 2));
                    
                    x_next = delta_U1 / fDelta;
                    y_next = delta_U2 / fDelta;
                    z_next = delta_U3 / fDelta;
                }

                // Случай с магнитным и постоянным электрическим полями
                if (mPlusEFields)
                {
                    double c1 = 2 * x - x_prev - alphaTemp * summ.bz * y_prev + alphaTemp * summ.by * z_prev;
                    double c2 = alphaTemp * summ.bx * x_prev + 2 * y - y_prev - alphaTemp * summ.bx * z_prev;
                    double c3 = - alphaTemp * summ.by * x_prev + alphaTemp * summ.bx * y_prev + 2 * z - z_prev;
                    
                    double fDelta = 1 + Math.Pow(alphaTemp, 2) * (Math.Pow(summ.bx, 2) + Math.Pow(summ.by, 2) + Math.Pow(summ.bz, 2));
                    
                    // Добавить + E_x(xi, yi, zi), + E_y(xi, yi, zi), + E_z(xi, yi, zi)
                    double delta_U1 = c1 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bx, 2)) + c2 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (alphaTemp * summ.bx * summ.bz - alphaTemp * Math.Pow(summ.by, 2));
                    double delta_U2 = - c1 * (alphaTemp * summ.bz + Math.Pow(alphaTemp, 2) * summ.bx * summ.by) + c2 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.by, 2)) + c3 * (alphaTemp * summ.bx + Math.Pow(alphaTemp, 2) * summ.by * summ.bz);
                    double delta_U3 = c1 * (Math.Pow(alphaTemp, 2) * summ.bx * summ.bz + alphaTemp * summ.by) - c2 * (alphaTemp * summ.bx - Math.Pow(alphaTemp, 2) * summ.by * summ.bz) + c3 * (1 + Math.Pow(alphaTemp, 2) * Math.Pow(summ.bz, 2));
                    
                    x_next = delta_U1 / fDelta;
                    y_next = delta_U2 / fDelta;
                    z_next = delta_U3 / fDelta;
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

                //Рисование
                //Vector3 pos = new Vector3((float)x * 100f, (float)y * 100f, (float)z * 100f);
                Vector3 pos = new Vector3((float)x, (float)y, (float)z);
                //Vector3 pos = new Vector3((float)x / 10f, (float)y / 10f, (float)z / 10f);
                AddLinePoint(ionPoints, 1, pos, ionsParams[ii].col);

                //Debug.Log("x diff = " + (x_next - x) + ", y diff = " + (y_next - y) + ", z diff = " + (z_next - z));

                // Переприсвоение перед новым циклом обсчёта
                x_prev = x; y_prev = y; z_prev = z;
                x = x_next; y = y_next; z = z_next;
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
    }
    
//     // Старый метод расчета точек движения заряженной частицы
//     private void CalculateIonTrajectory()
//     {
//         for (int ii = 0; ii < ionsParams.Count; ii++)
//         {
//             List<List<LinePoint>> ionPoints = new List<List<LinePoint>>();
//
//             double x_prev = p0.x / 1000f; double y_prev = p0.y / 1000f; double z_prev = p0.z / 1000f; //в метрах
//             double x = p1.x / 1000f; double y = p1.y / 1000f; double z = p1.z / 1000f; //в метрах
//
//             double h = Math.Pow(10, -6);
//             //double h = Math.Pow(10, -4);
//
//             double omega = 80 * Math.PI;
//             //double omega = 200 * Math.PI;
//
//             for (int j = 0; j < ionSteps; j++)
//             {
//                 // Подсчёт значения B
//                 DQuaternion M_i = new DQuaternion(x, y, z, 0, 0, 0, 0);
//
//                 DQuaternion summ = new DQuaternion();
//
//                 foreach (var coil in coilsP)
//                 {
//                     foreach (var ring in coil.rings)
//                     {
//                         R = coil.R;
//                         i = coil.i;
//
//                         ToLocalCartesian(ref M_i, ring);
//                         LocalCartesian2LocalCylindrical(M_i);
//                         CalculateB();
//                         M_i = LocalCylindrical2LocalCartesian(M_i);
//                         CalculateM(ref M_i);
//                         ToGlobalCartesian(ref M_i, ring);
//
//                         summ.x = M_i.x;
//                         summ.y = M_i.y;
//                         summ.z = M_i.z;
//
//                         summ.w += M_i.w; // * coil.turnsNum; // умножение на количество намоток
//                         summ.bx += M_i.bx;
//                         summ.by += M_i.by;
//                         summ.bz += M_i.bz;
//                     }
//                 }
//
//                 //double B_const = Math.Pow(10, -3);
//
//                 //Debug.Log(summ.w);
//
//                 //double C_ = ionsParams[ii].q * h * h / ionsParams[ii].m;
//
//                 // Вектор Bz должен быть сонаправлен оси z
//                 if (summ.w * summ.bz < 0)
//                     summ.w = summ.w * -1;
//
//                 double K_i_mag = (ionsParams[ii].q * h / (2 * ionsParams[ii].m)) * summ.w;
//                 double K_i_elmag = (ionsParams[ii].q * h / (2 * ionsParams[ii].m)) * summ.w * Math.Cos(omega * j * h);
//                 //double K_p = (C_ / 2 * h) * summ.w * Math.Cos(omega3); // Не хватает домножения на t
//
//
//                 // Расчет электрической составляющей
//                 //double E = 1; // В/м
//                 double E = 1d / 1000; // В/м
//                 //double L_i = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E * Math.Cos(omega * j * h);
//                 double L_i_mag = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E;
//                 double L_i_elmag = (ionsParams[ii].q * Math.Pow(h, 2) / ionsParams[ii].m) * E * Math.Cos(omega * j * h);
//
//                 //double C1 = 2 * x - x_prev - K_i * y_prev;
//                 //double C2 = 2 * y - y_prev + K_i * x_prev;
//                 //double delta_x = (C1 + K_i * C2);
//                 //double delta_y = (C2 - K_i * C1);
//                 //double fDelta = (1 + Math.Pow(K_i, 2));
//
//                 //Вычислим углы математически
//                 double cosGamma = Math.Abs(x * 0d + y * 0d + z * 1d) / (Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2)) * Math.Sqrt(Math.Pow(0d, 2) + Math.Pow(0d, 2) + Math.Pow(1d, 2))); // с вектором Vector3(0, 0, 1)
//                 double sinGamma = Math.Sqrt(1d - Math.Pow(cosGamma, 2));
//                 double cosBeta = Math.Abs(x * 1d + y * 0d + z * 0d) / (Math.Sqrt(Math.Pow(x, 2) + Math.Pow(y, 2) + Math.Pow(z, 2)) * Math.Sqrt(Math.Pow(1d, 2) + Math.Pow(0d, 2) + Math.Pow(0d, 2))); // с вектором Vector3(1, 0, 0)
//                 double sinBeta = Math.Sqrt(1d - Math.Pow(cosBeta, 2));
//
//                 double angGamma = Math.Acos(cosGamma) * (180d / Math.PI);
//                 double angBeta = Math.Acos(cosBeta) * (180d / Math.PI);
//
//                 Vector3 pointVec = new Vector3((float)x, (float)y, (float)z);
//                 //gammaAngle = Vector3.Angle(new Vector3(0, 0, 1), pointVec);
//                 //betaAngle = Vector3.Angle(new Vector3(1, 0, 0), pointVec);
//                 gammaAngle = Vector3.Angle(new Vector3(0, 0.66f, 0.33f), pointVec);
//                 betaAngle = Vector3.Angle(new Vector3(0.45f, 0.45f, 0), pointVec);
//
//                 double x_next = 0;
//                 double y_next = 0;
//                 double z_next = 0;
//
//                 // Слуйчас с магнитным полем, без электрического
//                 if (mField)
//                 {
//                     x_next = (2 * x - x_prev - K_i_mag * y_prev + K_i_mag * (2 * y - y_prev + K_i_mag * x_prev)) / (1 + Math.Pow(K_i_mag, 2));
//                     y_next = (2 * y - y_prev + K_i_mag * x_prev - K_i_mag * (2 * x - x_prev - K_i_mag * y_prev)) / (1 + Math.Pow(K_i_mag, 2));
//                     z_next = 2 * z - z_prev;
//                 }
//
//                 // Слуйчас с магнитным и постоянным электрическим полями
//                 if (mPlusEFields)
//                 {
//                     double fDelta = 1 + Math.Pow(K_i_mag, 2);
//                     double delta_x = 2 * x - x_prev - K_i_mag * y_prev + L_i_mag * sinGamma * cosBeta + K_i_mag * (2 * y - y_prev + K_i_mag * x_prev + L_i_mag * sinGamma * sinBeta);
//                     double delta_y = 2 * y - y_prev + K_i_mag * x_prev + L_i_mag * sinGamma * sinBeta - K_i_mag * (2 * x - x_prev - K_i_mag * y_prev + L_i_mag * sinGamma * cosBeta);
//                     x_next = delta_x / fDelta;
//                     y_next = delta_y / fDelta;
//                     z_next = 2 * z - z_prev + L_i_mag * cosGamma;
//                 }
//
//                 // Слуйчас с магнитным и переменным электрическим полями
//                 if (mPlusDEFields)
//                 {
//                     double fDelta = 1 + Math.Pow(K_i_elmag, 2);
//                     double delta_x = 2 * x - x_prev - K_i_elmag * y_prev + L_i_elmag * sinGamma * cosBeta + K_i_elmag * (2 * y - y_prev + K_i_elmag * x_prev + L_i_elmag * sinGamma * sinBeta);
//                     double delta_y = 2 * y - y_prev + K_i_elmag * x_prev + L_i_elmag * sinGamma * sinBeta - K_i_elmag * (2 * x - x_prev - K_i_elmag * y_prev + L_i_elmag * sinGamma * cosBeta);
//                     x_next = delta_x / fDelta;
//                     y_next = delta_y / fDelta;
//                     z_next = 2 * z - z_prev + L_i_elmag * cosGamma;
//                 }
//
//                 // Отладка и проверка векторов
//                 double x_dot = (x_next - x_prev) / (2 * h);
//                 double y_dot = (y_next - y_prev) / (2 * h);
//                 double z_dot = (z_next - z_prev) / (2 * h);
// //                Debug.Log("(x, y, z) = (" + x_next + ", " + y_next + ", " + z_next + ")" + ", (x_d, y_d, z_d) = (" + x_dot + ", " + y_dot + ", " + z_dot + ")");
//
//                 // Вчисление расстояния до самой удаленной точки для масштабирования
//                 double maxOfThree = Math.Max(x, Math.Max(y, z));
//                 if (maxOfThree > maxIonRange)
//                     maxIonRange = (float)maxOfThree;
//
//                 //Рисование
//                 //Vector3 pos = new Vector3((float)x * 100f, (float)y * 100f, (float)z * 100f);
//                 Vector3 pos = new Vector3((float)x, (float)y, (float)z);
//                 //Vector3 pos = new Vector3((float)x / 10f, (float)y / 10f, (float)z / 10f);
//                 AddLinePoint(ionPoints, 1, pos, ionsParams[ii].col);
//
//                 //Debug.Log("x diff = " + (x_next - x) + ", y diff = " + (y_next - y) + ", z diff = " + (z_next - z));
//
//                 // Переприсвоение перед новым циклом обсчёта
//                 x_prev = x; y_prev = y; z_prev = z;
//                 x = x_next; y = y_next; z = z_next;
//             }
//
//             Material _material;
//
//             switch (ii)
//             {
//                 case 1:
//                     _material = matLine3;
//                     break;
//                 case 2:
//                     _material = matLine2;
//                     break;
//                 default:
//                     _material = matLine1;
//                     break;
//             }
//
//             Debug.Log(ionPoints.Count);
//             CreateLinesMesh(ionPoints, _material);
//
//             //CreateCompressedLinesMesh(ionPoints, ref meshHolderIon);
//         }
//     }

    // Расчет визуализации магнитного поля
    private void CalculatePowerLine()
    {
        List<List<LinePoint>> points = new List<List<LinePoint>>();

        int gridStep = 0;

        var rings = new List<GameObject>();

        foreach (var coil in coilsP)
        {
            Grow(true, ref gridStep, points, coil);
            Grow(false, ref gridStep, points, coil);   
        }

        Material _material = matMag;

        Debug.Log(points.Count);
        CreateLinesMesh(points, _material);
    }

    private void Grow(bool dir, ref int gridStep, List<List<LinePoint>> points, Coil growCoil)
    {
        double localMeshStep = meshStep * 1000f;

        maxM = -1000;
        minM = 1000;

        //for (float xi = -5f; xi <= 5f; xi += (float) localMeshStep / 20f) // 30, 2.5
        for (float xi = -5f; xi <= 5f; xi += (float) localMeshStep / 20f) // 30, 2.5
        {
            //for (float yi = -5f; yi <= 5f; yi += (float) localMeshStep / 20f, gridStep++)
            for (float yi = -5f; yi <= 5f; yi += (float)localMeshStep / 20f, gridStep++)
            {
                // Вырезание окрестностей ноля
                //if ((xi > -4.15f) && ((xi < -3.85f)))
                //{
                //    xi = -xi;
                //    Debug.Log(xi);
                //}
                //if ((yi > -4.15f) && ((yi < -3.85f))) yi = -yi;

                float x_i = xi / 100f;
                float y_i = yi / 100f;
                float z_i = 0;

                var startPos = new Vector3(x_i, y_i, z_i);
                startPos = growCoil.transform.TransformDirection(startPos)+growCoil.transform.position;
                x_i = startPos.x;
                y_i = startPos.y;
                z_i = startPos.z;

                //float delta = 0.05f;

                float stepNumber = 0;
                float h = 0.1f;

                //Color lowColor = Color.blue;
                //Color highColor = Color.red;
                Color lowColor = new Color(204f / 255f, 242f / 255f, 255f / 255f);
                Color highColor = new Color(1f, 153f / 255f, 184f / 255f);
                lowColor.a = colorIntensity;
                highColor.a = colorIntensity;

                while ((x_i >= (-xFieldSize / 2f) * 10f) && (x_i <= (xFieldSize / 2f) * 10f) 
                    && (y_i >= (-yFieldSize / 2f) * 10f) && (y_i <= (yFieldSize / 2f) * 10f) 
                    && (z_i >= (-zFieldSize / 2f) * 10f) && (z_i <= (zFieldSize / 2f) * 10f) 
                    && (stepNumber < maxSteps))
                {
                    stepNumber++;

                    // Подсчёт значения B
                    DQuaternion M_i = new DQuaternion(x_i, y_i, z_i, 0, 0, 0, 0);

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

                    // ПЕРЕПИСАТЬ — ЭТО НЕ СОВСЕМ ПРАВДА
                    if ((summ.w > maxM) && (!double.IsInfinity(summ.w))) maxM = summ.w;
                    if (summ.w < minM) minM = summ.w;

                    double localMaxM = maxM - minM;

                    // Визуализация шага всторону
                    var pos = new Vector3((float) summ.x, (float) summ.y, (float) summ.z);
                    var color = Color.Lerp(lowColor, highColor, (float) (10 * summ.w / localMaxM));
                    AddLinePoint(points, gridStep, pos, color);

                    //SpawnGoSegment(line_i_right, pos, summ, color);
                    
                    // Расчет точек смещения
                    double alpha_i = Math.Sqrt(Math.Pow(summ.bx, 2) + Math.Pow(summ.by, 2) + Math.Pow(summ.bz, 2));
                    x_i = x_i + ((h * (float) summ.bx) / (float) alpha_i)*(dir?1:-1);
                    y_i = y_i + ((h * (float) summ.by) / (float) alpha_i)*(dir?1:-1);
                    z_i = z_i + ((h * (float) summ.bz) / (float) alpha_i)*(dir?1:-1);
                }
            }
        }
    }

    private void SpawnGoSegment(Vector3 pos, DQuaternion summ, Color color)
    {
        GameObject line_i_right = Instantiate(lineSource);

        line_i_right.transform.parent = fieldSource.transform;
        line_i_right.transform.position = pos;
        Vector3 orientationPointRight =
            new Vector3((float) summ.bx * 10000, (float) summ.@by * 10000, (float) summ.bz * 10000);
        line_i_right.transform.LookAt(orientationPointRight.normalized + line_i_right.transform.position);
        Renderer rnd2 = line_i_right.GetComponentInChildren<Renderer>();
        MaterialPropertyBlock props2 = new MaterialPropertyBlock();
        //props2.SetColor("_Color", color);
        props2.SetColor("_Color", color);
        rnd2.SetPropertyBlock(props2);
        rnd2.gameObject.layer = 8;
    }

    //[ContextMenu("TestLineDraw")]
    //public void TestLineDraw()
    //{
    //    var points = new List<List<LinePoint>>();


    //    for (int x = 0; x < 5; x++)
    //    {
    //        points.Add(new List<LinePoint>());

    //        for (int y = 0; y < 30; y++)
    //        {
    //            points[x].Add(new LinePoint()
    //            {
    //                color = Color.Lerp(Color.red, Color.green, Random.Range(0f, 1f)),
    //                pos = new Vector3(Mathf.Sin(x), Mathf.Sin(y), y)
    //            });
    //        }
    //    }


    //    CreateLinesMesh(points);
    //}

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


    //private void CreateLinesMesh(List<List<LinePoint>> points, ref GameObject meshHolderObj)
    //{
    //    for (int i = points.Count - 1; i >= 0; i--)
    //    {
    //        if (points[i].Count < 2)
    //        {
    //            points.RemoveAt(i);
    //            continue;
    //        }
    //    }

    //    List<Vector3> positions = new List<Vector3>();
    //    List<int> indexes = new List<int>();
    //    List<Color> colors = new List<Color>();

    //    var indexCounter = 0;
    //    for (int x = 0; x < points.Count; x++)
    //    {
    //        for (int y = 1; y < points[x].Count; y++)
    //        {
    //            positions.Add(points[x][y - 1].pos);
    //            positions.Add(points[x][y].pos);
    //            colors.Add(points[x][y - 1].color);
    //            colors.Add(points[x][y].color);
    //            indexes.Add(indexCounter++);
    //            indexes.Add(indexCounter++);
    //        }
    //    }

    //    var mesh = new Mesh();
    //    mesh.bounds = new Bounds(Vector3.zero, Vector3.one * 100500);
    //    mesh.indexFormat = IndexFormat.UInt32;
    //    mesh.vertices = positions.ToArray();
    //    mesh.colors = colors.ToArray();
    //    mesh.SetIndices(indexes.ToArray(), MeshTopology.Lines, 0, false);

    //    meshHolderObj = new GameObject("MeshHolder");
    //    var filter = meshHolderObj.AddComponent<MeshFilter>();
    //    filter.sharedMesh = mesh;
    //    var renderer = meshHolderObj.AddComponent<MeshRenderer>();
    //    //var shader = Shader.Find("Particles/Standard Unlit");
    //    //var mat = new Material(shader);
    //    mat.color = Color.gray;
    //    renderer.sharedMaterial = mat;
    //    renderer.receiveShadows = false;
    //    renderer.shadowCastingMode = ShadowCastingMode.Off;
    //}

    private void CreateCompressedLinesMesh(List<List<LinePoint>> points, ref GameObject meshHolderObj)
    {
        for (int i = points.Count - 1; i >= 0; i--)
        {
            if (points[i].Count < 2)
            {
                points.RemoveAt(i);
                continue;
            }
        }

        List<Vector3> positions = new List<Vector3>();
        List<int> indexes = new List<int>();
        List<Color> colors = new List<Color>();

        var indexCounter = 0;

        float scaleCoeff = 1f;

        //if ((maxIonRange * 100f) > (xFieldSize / 2f))
        //    scaleCoeff = (xFieldSize / 2f) / (maxIonRange * 100f);

        for (int x = 0; x < points.Count; x++)
        {
            for (int y = 1; y < points[x].Count; y++)
            {

                Vector3 tempVecPos1 = new Vector3(points[x][y - 1].pos.x * scaleCoeff, points[x][y - 1].pos.y * scaleCoeff, points[x][y - 1].pos.z * scaleCoeff);
                Vector3 tempVecPos2 = new Vector3(points[x][y].pos.x * scaleCoeff, points[x][y].pos.y * scaleCoeff, points[x][y].pos.z * scaleCoeff);

                positions.Add(tempVecPos1);
                positions.Add(tempVecPos2);
                colors.Add(points[x][y - 1].color);
                colors.Add(points[x][y].color);
                indexes.Add(indexCounter++);
                indexes.Add(indexCounter++);
            }
        }

        var mesh = new Mesh();
        mesh.bounds = new Bounds(Vector3.zero, Vector3.one * 100500);
        mesh.indexFormat = IndexFormat.UInt32;
        mesh.vertices = positions.ToArray();
        mesh.colors = colors.ToArray();
        mesh.SetIndices(indexes.ToArray(), MeshTopology.Lines, 0, false);

        meshHolderObj = new GameObject("MeshHolder");
        var filter = meshHolderObj.AddComponent<MeshFilter>();
        filter.sharedMesh = mesh;
        var renderer = meshHolderObj.AddComponent<MeshRenderer>();
        var shader = Shader.Find("Particles/Standard Unlit");
        var mat = new Material(shader);
        mat.color = Color.gray;
        renderer.sharedMaterial = mat;
        renderer.receiveShadows = false;
        renderer.shadowCastingMode = ShadowCastingMode.Off;
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
}
