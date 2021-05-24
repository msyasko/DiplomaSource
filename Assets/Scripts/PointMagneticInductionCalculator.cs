using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Threading;
using UnityEngine;
using Debug = UnityEngine.Debug;

[Serializable]
public struct DQuaternion
{
    public double x;
    public double y;
    public double z;
    public double w;
    public double bx;
    public double by;
    public double bz;

    public DQuaternion(double _x, double _y, double _z, double _w, double _bx, double _by, double _bz) 
    {
        x = _x;
        y = _y;
        z = _z;
        w = _w;
        bx = _bx;
        by = _by;
        bz = _bz;
    }
}

public struct DVector3
{
    public double x;
    public double y;
    public double z;

    public DVector3(double _x, double _y, double _z)
    {
        x = _x;
        y = _y;
        z = _z;
    }
}

[Serializable]
public class DQuaternionList
{
    private List<DQuaternion> l = new List<DQuaternion>();
    private Dictionary<int, DQuaternion> _dic = new Dictionary<int, DQuaternion>(); 

    public int currentWorkersCount => _currentWorkersCount;

    public int Count
    {
        get { return l.Count; }
    }

    public int CountDictionary
    {
        get
        {
            int countL = 0;
            lock (_dic)
            {
                countL = _dic.Count;
            }
            return countL;
        }
    }

    private int _currentWorkersCount;
    private object _lockObject = new object();

    public void AddSafely(DQuaternion value)
    {
        lock (l)
        {
            l.Add(value);
        }
    }

    public void AddSafely(DQuaternion value, int index)
    {
        lock (_dic)
        {
            _dic[index] = value;
        }
    }

    public void RegisterWorker()
    {
        lock (_lockObject)
        {
            _currentWorkersCount++;
        }
    }

    public void UnRegisterWorker()
    {
        lock (_lockObject)
        {
            _currentWorkersCount--;
        }
    }

    internal DQuaternion GetPoint(int j)
    {
        lock (l)
        {
            return l[j];
        }
    }

    internal DQuaternion GetPointByIndex(int index)
    {
        lock (_dic)
        {
            DQuaternion i;
            if (_dic.TryGetValue(index, out i))
            return i;
            Debug.LogErrorFormat("Wrong index : {0}", index);
            return new DQuaternion(); //_dic[index];
        }
    }
}

public class PointMagneticInductionCalculator : MonoBehaviour
{

    public ParticleSystem fieldPoints;
    public List<ParticleSystem.Particle> particleList = new List<ParticleSystem.Particle>();

    private DVector3 B = new DVector3();
    public DQuaternionList MList = new DQuaternionList();

    private double z; //
    public double R = 1.0f; // радиус контура S
    public double i = 1.0f; // величина силы тока, протекающего по контуру
    public double O = 0.0f; // центр контура S
    private double ro; // расстояние от точки M до оси Oz
    public double dl = 0.0f; // элементарный участок контура
    public double doubleAlpha = 30.0f; // угол между двумя элементарными участками контура, расположенными симметрично
    private double beta;
    private double r; // расстояние от точки M до элементарных участков контура dl
    public double mu0; // магнитная постоянная

    private double K;
    private double L;
    private double B_ro;
    private double B_z;

    public double maxM;
    public double minM;

    //public double meshRange = 2.0f; // размер сетки
    public double meshStep = 20f; // шаг обсчёта

    public GameObject lineSource;
    public GameObject fieldSource;
    public GameObject pipeSource;
    private GameObject[] lineArray;
    private Vector3[] magneticFieldArray;

    // Задание параметров рассматриваемой области
    public float xFieldSize;
    public float yFieldSize;
    public float zFieldSize;
    public float inductorRadius;
    public float windingNum;

    public delegate double integratedFunction(double x); //Объявление делегата интегрируемой ф-ции

    // Use this for initialization
    void Start()
    {
        Stopwatch sw = new Stopwatch();
        sw.Start();

        r = (Math.Pow(R, 2) + Math.Pow(ro, 2) - 2 * ro * R * Math.Cos(doubleAlpha / 2) + Math.Pow(z, 2));
        mu0 = (4 * Math.PI * Math.Pow(10, -7));

        beta = ((doubleAlpha / 2 - Math.PI) / 2);

        int workerThreads, completionPortThreads;
        ThreadPool.GetMaxThreads(out workerThreads, out completionPortThreads);
        Debug.LogWarningFormat("Доступно потоков : {0}", workerThreads);
        //new Action(CalculateMagneticField).BeginInvoke(null, null);
        CalculateMagneticField();

        // Отрисовка
        //InstantiateVisualisationObjects();
        //StartCoroutine(SetupVisualistaionSystem());

        CalculatePowerLine();

        //StartCoroutine(SetupParticleSystem());
        sw.Stop();
        Debug.LogWarningFormat("Time to calculate : {0}", sw.ElapsedMilliseconds);
    }

    // Update is called once per frame
    //void Update()
    //{
    //}

    //=============================================================
    // Подготовка рабочей области и заполнение массива точек
    //=============================================================
    void instantiatePointField()
    {
        for (float x = -xFieldSize / 2f; x <= xFieldSize / 2f; x += (float)meshStep)
            for (float y = -yFieldSize / 2f; y <= yFieldSize / 2f; y += (float)meshStep)
                for (float z = -zFieldSize / 2f; z <= zFieldSize / 2f; z += (float)meshStep)
                {
                    DQuaternion M = new DQuaternion(x, y, z, 0, 0, 0, 0);
                    MList.AddSafely(M);
                }
    }

    //======================================================
    // Подготовка необходимых заранее вычислений
    //======================================================
    #region 
    // Приведение точки к локальным сферическим коордтинатам
    double transformPointToLocal(double beta)
    {
        return (1 / Math.Sqrt(1 - kSquared(ro, R, z) * Math.Pow(Math.Sin(beta), 2)));
    }

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

    // Вычисление значения интеграла K
    void CalculateK()
    {
        K = integration(integratedFunctionK, 0, Math.PI / 2);
    }

    // Вычисление значения интеграла L
    void CalculateL()
    {
        L = integration(integratedFunctionL, 0, Math.PI / 2);
    }

    // Вычисление векторного потенциала A в точке M
    double CalculateA()
    {
        return ((mu0 * i / (2 * Math.PI)) * (Math.Sqrt(R / ro)) * ((2 - kSquared(ro, R, z) * K / Math.Sqrt(kSquared(ro, R, z)) - (2 * L) / Math.Sqrt(kSquared(ro, R, z)))));
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
    void InstantiateVisualisationObjects()
    {
        int lineCounter = (int)(xFieldSize/meshStep + 1) * (int)(yFieldSize / meshStep + 1) * (int)(zFieldSize / meshStep + 1);
        lineArray = new GameObject[lineCounter];

        for (int i = 0; i < lineCounter - 1; i++)
        {
            lineArray[i] = Instantiate(lineSource);
            lineArray[i].transform.parent = fieldSource.transform;
        }
    }

    //======================================================
    // Проведение основных вычислений потенциала в точке M
    //======================================================
    #region
    // Переход из локальной ДСК в локальную ЦСК
    void TranformCaCStoCiCS(ref DQuaternion M)
    {
        z = M.z;
        ro = Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2));
    }

    // Переход из локальной ДСК в локальную ЦСК, если катушка не в 000
    void TranformCaCStoCiCS(ref DQuaternion M, Vector3 v3)
    {
        z = Math.Abs(M.z - v3.z);
        ro = Math.Sqrt(Math.Pow(Math.Abs(M.x - v3.x), 2) + Math.Pow(Math.Abs(M.y - v3.y), 2));
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
    void TranformCiCStoCaCS(ref DQuaternion M)
    {
        B = new DVector3(B_ro * M.x / (Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2))), B_ro * M.y / (Math.Sqrt(Math.Pow(M.x, 2) + Math.Pow(M.y, 2))), B_z);

        M.bx = B.x;
        M.by = B.y;
        M.bz = B.z;
    }

    // Перегрузка перехода из локальной ЦСК в локальную ДСК
    DQuaternion TranformCiCStoCaCS(DQuaternion dQM)
    {
        B = new DVector3(B_ro * dQM.x / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))), B_ro * dQM.y / (Math.Sqrt(Math.Pow(dQM.x, 2) + Math.Pow(dQM.y, 2))), B_z);
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

        // Запускается обход всех M
        for (double xi = -xFieldSize/2f; xi <= xFieldSize/2f; xi += meshStep)
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
        //int i = index;
        for (double yi = -yFieldSize / 2f; yi <= yFieldSize / 2f; yi += meshStep)
        {
            for (double zi = -zFieldSize / 2f; zi <= zFieldSize / 2f; zi += meshStep)
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

        TranformCaCStoCiCS(ref M);

        CalculateB();

        M = TranformCiCStoCaCS(M);

        CalculateM(ref M);

        if (!double.IsNaN(M.w) || !double.IsInfinity(M.w))
        {
            // Обрубаются близкие к бесконечности значения
            if (M.w < 1000)
            {
                MList.AddSafely(M);
                //MList.AddSafely(M, index);
                if ((M.w > maxM) && (!double.IsInfinity(M.w))) maxM = M.w;
                if (M.w < minM) minM = M.w;
            }
        }
        else
        {
            Debug.Log("( " + M.x + " , " + M.y + " , " + M.z + " , " + M.w);
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

        Color lowColor = Color.blue;
        lowColor.a = .1f;

        MaterialPropertyBlock props = new MaterialPropertyBlock();

        for (int j = 0; j < counter - 1; j++)
        {
            double localMaxM = maxM - minM;
            var tempM = MList.GetPoint(j); //MList.GetPointByIndex(j);

            lineArray[j].transform.position = new Vector3((float)tempM.x, (float)tempM.y, (float)tempM.z);
            Vector3 orientationPoint = new Vector3((float)tempM.bx * 10000, (float)tempM.by * 10000, (float)tempM.bz * 10000); 
            lineArray[j].transform.LookAt(orientationPoint.normalized + lineArray[j].transform.position);

            Renderer rnd = lineArray[j].GetComponentInChildren<Renderer>();
            props.SetColor("_Color", Color.Lerp(lowColor, Color.red, (float)(10 * tempM.w / localMaxM)));
            rnd.SetPropertyBlock(props);
            rnd.gameObject.layer = 8;

            //Debug.Log("Color: " + prt.color + ", Value: " + (tempM.w / localMaxM));

        }

        //fieldPoints.SetParticles(particleList.ToArray(), particleList.Count);
        sw.Stop();
        Debug.LogWarningFormat("Threaded work time seems to : {0}", sw.ElapsedMilliseconds);
    
    
    }

    // Пробы расчета и визуализации электромагнитных линий методом приближения
    //private void CalculatePowerLine()
    //{
    //    float x_i;
    //    float y_i;
    //    float z_i;

    //    float h = (float) meshStep;

    //    int counter = MList.Count;

    //    Color lowColor = Color.blue;
    //    lowColor.a = .1f;

    //    MaterialPropertyBlock props = new MaterialPropertyBlock();

    //    for (int j = 0; j < counter - 1; j++)
    //    {
    //        x_i = (float) MList.GetPoint(j).x;
    //        y_i = (float) MList.GetPoint(j).y;
    //        z_i = (float) MList.GetPoint(j).z;

    //        double localMaxM = maxM - minM;
    //        var tempM = MList.GetPoint(j); //MList.GetPointByIndex(j);

    //        GameObject line_i = Instantiate(lineSource);
    //        line_i.transform.parent = pipeSource.transform;

    //        double alpha_i = Math.Sqrt(Math.Pow(MList.GetPoint(j).bx, 2) + Math.Pow(MList.GetPoint(j).by, 2) + Math.Pow(MList.GetPoint(j).bz, 2));
    //        x_i = x_i + (h * (float)MList.GetPoint(j).bx) / (float)alpha_i;
    //        y_i = y_i + (h * (float)MList.GetPoint(j).by) / (float)alpha_i;
    //        z_i = z_i + (h * (float)MList.GetPoint(j).bz) / (float)alpha_i;

    //        line_i.transform.position = new Vector3((float)MList.GetPoint(j).x, (float)MList.GetPoint(j).y, (float)MList.GetPoint(j).z);
    //        Vector3 orientationPoint = new Vector3(x_i, y_i, z_i);
    //        line_i.transform.LookAt(orientationPoint.normalized + line_i.transform.position);

    //        Renderer rnd = line_i.GetComponentInChildren<Renderer>();
    //        props.SetColor("_Color", Color.Lerp(lowColor, Color.red, (float)(10 * tempM.w / localMaxM)));
    //        rnd.SetPropertyBlock(props);
    //        rnd.gameObject.layer = 8;

    //        //Debug.Log("Color: " + prt.color + ", Value: " + (tempM.w / localMaxM));

    //    }

    //    //for (float zi = 0; zi <= zFieldSize / 2f; zi += (float) meshStep)
    //    //{
    //        //for (float yi = 0; yi <= yFieldSize / 2f; yi += (float) meshStep)
    //        //{
    //        //    float x_i = 0.0f;
    //        //    float y_i = yi / 100f;
    //        //    float z_i = zi / 100f;

    //        //    float h = 0.1f;

    //        //    Color lowColor = Color.blue;
    //        //    lowColor.a = .1f;

    //        //    while ((x_i < xFieldSize / 100f) && (y_i < yFieldSize / 100f) && (z_i < zFieldSize / 100f) && (z_i >= 0))
    //        //    {
    //        //        // Подсчёт значения B
    //        //        DQuaternion M_i = new DQuaternion(x_i, y_i, z_i, 0, 0, 0, 0);

    //        //        TranformCaCStoCiCS(ref M_i);
    //        //        CalculateB();
    //        //        M_i = TranformCiCStoCaCS(M_i);
    //        //        CalculateM(ref M_i);

    //        //        // Визуализация
    //        //        GameObject line_i = Instantiate(lineSource);

    //        //        line_i.transform.parent = pipeSource.transform;

    //        //        line_i.transform.position = new Vector3((float)M_i.x, (float)M_i.y, (float)M_i.z);
    //        //        Vector3 orientationPoint = new Vector3((float)M_i.bx * 10000, (float)M_i.by * 10000, (float)M_i.bz * 10000);
    //        //        line_i.transform.LookAt(orientationPoint.normalized + line_i.transform.position);

    //        //        double localMaxM = maxM - minM;

    //        //        Renderer rnd = line_i.GetComponentInChildren<Renderer>();
    //        //        MaterialPropertyBlock props = new MaterialPropertyBlock();
    //        //        props.SetColor("_Color", Color.Lerp(lowColor, Color.red, (float)(10 * M_i.w / localMaxM)));
    //        //        rnd.SetPropertyBlock(props);
    //        //        rnd.gameObject.layer = 8;

    //        //        double alpha_i = Math.Sqrt(Math.Pow(M_i.bx, 2) + Math.Pow(M_i.by, 2) + Math.Pow(M_i.bz, 2));
    //        //        x_i = x_i + (h * (float)M_i.bx) / (float)alpha_i;
    //        //        y_i = y_i + (h * (float)M_i.by) / (float)alpha_i;
    //        //        z_i = z_i + (h * (float)M_i.bz) / (float)alpha_i;

    //        //        // Дублирование
    //        //        ///////
    //        //        //GameObject line_i_1 = Instantiate(line_i);
    //        //        //line_i_1.transform.parent = pipeSource.transform;
    //        //        //line_i_1.transform.position = new Vector3((float)-M_i.x, (float)M_i.y, (float)M_i.z);
    //        //        //orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)M_i.by * 10000, (float)M_i.bz * 10000);
    //        //        //line_i_1.transform.LookAt(orientationPoint.normalized + line_i_1.transform.position);
    //        //        //rnd = line_i_1.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_2 = Instantiate(line_i);
    //        //        //line_i_2.transform.parent = pipeSource.transform;
    //        //        //line_i_2.transform.position = new Vector3((float)-M_i.x, (float)-M_i.y, (float)M_i.z);
    //        //        //orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)-M_i.by * 10000, (float)M_i.bz * 10000);
    //        //        //line_i_2.transform.LookAt(orientationPoint.normalized + line_i_2.transform.position);
    //        //        //rnd = line_i_2.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_3 = Instantiate(line_i);
    //        //        //line_i_3.transform.parent = pipeSource.transform;
    //        //        //line_i_3.transform.position = new Vector3((float)M_i.x, (float)-M_i.y, (float)M_i.z);
    //        //        //orientationPoint = new Vector3((float)M_i.bx * 10000, (float)-M_i.by * 10000, (float)M_i.bz * 10000);
    //        //        //line_i_3.transform.LookAt(orientationPoint.normalized + line_i_3.transform.position);
    //        //        //rnd = line_i_3.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_4 = Instantiate(line_i);
    //        //        //line_i_4.transform.parent = pipeSource.transform;
    //        //        //line_i_4.transform.position = new Vector3((float)M_i.x, (float)M_i.y, (float)-M_i.z);
    //        //        //orientationPoint = new Vector3((float)M_i.bx * 10000, (float)M_i.by * 10000, (float)-M_i.bz * 10000);
    //        //        //line_i_4.transform.LookAt(orientationPoint.normalized + line_i_4.transform.position);
    //        //        //rnd = line_i_4.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_5 = Instantiate(line_i);
    //        //        //line_i_5.transform.parent = pipeSource.transform;
    //        //        //line_i_5.transform.position = new Vector3((float)-M_i.x, (float)M_i.y, (float)-M_i.z);
    //        //        //orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)M_i.by * 10000, (float)-M_i.bz * 10000);
    //        //        //line_i_5.transform.LookAt(orientationPoint.normalized + line_i_5.transform.position);
    //        //        //rnd = line_i_5.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_6 = Instantiate(line_i);
    //        //        //line_i_6.transform.parent = pipeSource.transform;
    //        //        //line_i_6.transform.position = new Vector3((float)-M_i.x, (float)-M_i.y, (float)-M_i.z);
    //        //        //orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)-M_i.by * 10000, (float)-M_i.bz * 10000);
    //        //        //line_i_6.transform.LookAt(orientationPoint.normalized + line_i_6.transform.position);
    //        //        //rnd = line_i_6.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;

    //        //        //GameObject line_i_7 = Instantiate(line_i);
    //        //        //line_i_7.transform.parent = pipeSource.transform;
    //        //        //line_i_7.transform.position = new Vector3((float)M_i.x, (float)-M_i.y, (float)-M_i.z);
    //        //        //orientationPoint = new Vector3((float)M_i.bx * 10000, (float)-M_i.by * 10000, (float)-M_i.bz * 10000);
    //        //        //line_i_7.transform.LookAt(orientationPoint.normalized + line_i_7.transform.position);
    //        //        //rnd = line_i_7.GetComponentInChildren<Renderer>();
    //        //        //rnd.SetPropertyBlock(props);
    //        //        //rnd.gameObject.layer = 8;
    //        //        /////////
    //        //    }
    //        //}
    //    //}
    //}

    private void CalculatePowerLine()
    {

        // Поменял xi и zi

        for (float xi = 0; xi <= xFieldSize / 2f; xi += (float)meshStep)
        {
            for (float yi = 0; yi <= yFieldSize / 2f; yi += (float)meshStep)
            {
                float x_i = xi / 100f;
                float y_i = yi / 100f;
                float z_i = 0.0f;

                float h = 0.1f;

                Color lowColor = Color.blue;
                lowColor.a = .1f;

                while ((x_i < xFieldSize / 100f) && (y_i < yFieldSize / 100f) && (z_i < zFieldSize / 100f) && (z_i >= 0))
                {
                    // Подсчёт значения B
                    DQuaternion M_i = new DQuaternion(x_i, y_i, z_i, 0, 0, 0, 0);

                    TranformCaCStoCiCS(ref M_i, pipeSource.gameObject.transform.position);
                    CalculateB();
                    M_i = TranformCiCStoCaCS(M_i);
                    CalculateM(ref M_i);

                    // Визуализация
                    GameObject line_i = Instantiate(lineSource);

                    line_i.transform.parent = fieldSource.transform;

                    line_i.transform.position = new Vector3((float)M_i.x, (float)M_i.y, (float)M_i.z);
                    Vector3 orientationPoint = new Vector3((float)M_i.bx * 10000, (float)M_i.by * 10000, (float)M_i.bz * 10000);
                    line_i.transform.LookAt(orientationPoint.normalized + line_i.transform.position);

                    double localMaxM = maxM - minM;

                    Renderer rnd = line_i.GetComponentInChildren<Renderer>();
                    MaterialPropertyBlock props = new MaterialPropertyBlock();
                    props.SetColor("_Color", Color.Lerp(lowColor, Color.red, (float)(10 * M_i.w / localMaxM)));
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    double alpha_i = Math.Sqrt(Math.Pow(M_i.bx, 2) + Math.Pow(M_i.by, 2) + Math.Pow(M_i.bz, 2));
                    x_i = x_i + (h * (float)M_i.bx) / (float)alpha_i;
                    y_i = y_i + (h * (float)M_i.by) / (float)alpha_i;
                    z_i = z_i + (h * (float)M_i.bz) / (float)alpha_i;

                    // Дублирование
                    ///////
                    GameObject line_i_1 = Instantiate(line_i);
                    line_i_1.transform.parent = fieldSource.transform;
                    line_i_1.transform.position = new Vector3((float)-M_i.x, (float)M_i.y, (float)M_i.z);
                    orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)M_i.by * 10000, (float)M_i.bz * 10000);
                    line_i_1.transform.LookAt(orientationPoint.normalized + line_i_1.transform.position);
                    rnd = line_i_1.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_2 = Instantiate(line_i);
                    line_i_2.transform.parent = fieldSource.transform;
                    line_i_2.transform.position = new Vector3((float)-M_i.x, (float)-M_i.y, (float)M_i.z);
                    orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)-M_i.by * 10000, (float)M_i.bz * 10000);
                    line_i_2.transform.LookAt(orientationPoint.normalized + line_i_2.transform.position);
                    rnd = line_i_2.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_3 = Instantiate(line_i);
                    line_i_3.transform.parent = fieldSource.transform;
                    line_i_3.transform.position = new Vector3((float)M_i.x, (float)-M_i.y, (float)M_i.z);
                    orientationPoint = new Vector3((float)M_i.bx * 10000, (float)-M_i.by * 10000, (float)M_i.bz * 10000);
                    line_i_3.transform.LookAt(orientationPoint.normalized + line_i_3.transform.position);
                    rnd = line_i_3.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_4 = Instantiate(line_i);
                    line_i_4.transform.parent = fieldSource.transform;
                    line_i_4.transform.position = new Vector3((float)M_i.x, (float)M_i.y, (float)-M_i.z);
                    orientationPoint = new Vector3((float)M_i.bx * 10000, (float)M_i.by * 10000, (float)-M_i.bz * 10000);
                    line_i_4.transform.LookAt(orientationPoint.normalized + line_i_4.transform.position);
                    rnd = line_i_4.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_5 = Instantiate(line_i);
                    line_i_5.transform.parent = fieldSource.transform;
                    line_i_5.transform.position = new Vector3((float)-M_i.x, (float)M_i.y, (float)-M_i.z);
                    orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)M_i.by * 10000, (float)-M_i.bz * 10000);
                    line_i_5.transform.LookAt(orientationPoint.normalized + line_i_5.transform.position);
                    rnd = line_i_5.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_6 = Instantiate(line_i);
                    line_i_6.transform.parent = fieldSource.transform;
                    line_i_6.transform.position = new Vector3((float)-M_i.x, (float)-M_i.y, (float)-M_i.z);
                    orientationPoint = new Vector3((float)-M_i.bx * 10000, (float)-M_i.by * 10000, (float)-M_i.bz * 10000);
                    line_i_6.transform.LookAt(orientationPoint.normalized + line_i_6.transform.position);
                    rnd = line_i_6.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;

                    GameObject line_i_7 = Instantiate(line_i);
                    line_i_7.transform.parent = fieldSource.transform;
                    line_i_7.transform.position = new Vector3((float)M_i.x, (float)-M_i.y, (float)-M_i.z);
                    orientationPoint = new Vector3((float)M_i.bx * 10000, (float)-M_i.by * 10000, (float)-M_i.bz * 10000);
                    line_i_7.transform.LookAt(orientationPoint.normalized + line_i_7.transform.position);
                    rnd = line_i_7.GetComponentInChildren<Renderer>();
                    rnd.SetPropertyBlock(props);
                    rnd.gameObject.layer = 8;
                    /////////
                }
            }
        }
    }

    //// Настройка визуализации
    //IEnumerator SetupParticleSystem()
    //{
    //    Stopwatch sw = new Stopwatch();
    //    sw.Start();
    //    yield return null;
    //    //yield return new WaitForSeconds(60);
    //    while (MList.currentWorkersCount > 0)
    //    {
    //        yield return null;
    //    }

    //    //yield return new WaitForSeconds(1f);

    //    //while (MList.currentWorkersCount > 0)
    //    //{
    //    //    yield return null;
    //    //}

    //    int counter =  MList.Count;//MList.CountDictinary;

    //    ParticleSystem.Particle prt = new ParticleSystem.Particle();
    //    prt.size = 0.05f;
    //    Color lowColor = Color.blue;
    //    lowColor.a = .1f;

    //    for (int j = 0; j < counter - 1; j++)
    //    {
    //        double localMaxM = maxM - minM;

    //        var tempM = MList.GetPoint(j); //MList.GetPointByIndex(j);

    //        prt.position = new Vector3((float)tempM.x, (float)tempM.y, (float)tempM.z);
            
    //        prt.color = Color.Lerp(lowColor, Color.red, (float)(10 * tempM.w / localMaxM));

    //        //Debug.Log("Color: " + prt.color + ", Value: " + (tempM.w / localMaxM));

    //        particleList.Add(prt); 
    //    }

    //    fieldPoints.SetParticles(particleList.ToArray(), particleList.Count);
    //    sw.Stop();
    //    Debug.LogWarningFormat("Threaded work time seems to : {0}", sw.ElapsedMilliseconds);
    //}
    #endregion



    //private void btCalc_Click(object sender, EventArgs e)  //Вызов метода интегрирования
    //{
    //    integratedFunction fun = new integratedFunction(имя_интегрируемой_функции);
    //    MessageBox.Show(integring(fun, нижний_предел_интегрирования, верхний_предел).ToString());
    //}
}
