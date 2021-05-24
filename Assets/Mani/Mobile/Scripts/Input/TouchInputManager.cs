using System;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;
using UnityEngine.EventSystems;

namespace Gestures
{
    public class TouchInputManager : MonoBehaviour
    {
        public event Action<Vector2> ETap;
        public event Action<Vector2> EDoubleTap;

        private readonly TouchProtector[] touchProtectorsP = new TouchProtector[10];
        private readonly PrevTouchData[] prevTouchDatasP = new PrevTouchData[10];
        
        public List<ExtendedTouch> extendedTouches;

        private static TouchInputManager instanceS;
        private float referenceDpi = 320;
        public int cojointDistance = 15;
        public float deadZone = 0.1f;

        public float tapTimeout = 0.2f;
        public float doubleTapTimeout = 0.2f;
        public float tapMaxRadius = 0.05f;

        public Texture touchCircle;
        public LayerMask uiMask;

        public static TouchInputManager Instance { get { return instanceS; } }

        public static float Dpi
        {
            get
            {
#if UNITY_EDITOR
                return Instance.referenceDpi;
#else
                return Screen.dpi;
# endif
            }
        }

        void Awake()
        {
            if (Instance)
                return;

            instanceS = this;
        }


        void Update()
        {
            extendedTouches.RemoveAll(t => t.phase == TouchPhase.Ended);
            List<ExtendedTouch> newTouches = new List<ExtendedTouch>();

            List<Touch> touches = Input.touches.Where(t=>t.fingerId<9).ToList();
            RegisterPrevTouches(touches);
            ApplyDeadZones(touches);

            FindCojoined(touches, newTouches, extendedTouches);
            FindTwist(touches, newTouches, extendedTouches, false);
            FindPinch(touches, newTouches, extendedTouches);
            MarkAsUsed(newTouches, touches);
            FindSingle(touches, newTouches, extendedTouches);
            FindEnded(newTouches);

            #region Test

            //if (extendedTouches.Count > 0)
            //{
            //    Debug.Log(Time.frameCount);
            //    foreach (var t in extendedTouches)
            //    {
            //        Debug.Log(t.type);
            //        Debug.Log(t.phase);
            //        Debug.Log(t.GetHashCode());
            //    }
            //}

            //if (extendedTouches.Count>0)
            //{
            //    int cj = 0;
            //    int cM = 0;
            //    int p = 0;
            //    int s = 0;
            //    int sM = 0;

            //    Debug.Log(Time.frameCount);
            //    foreach (var t in extendedTouches)
            //    {
            //        switch (t.type)
            //        {
            //            case TouchType.Conjoined:
            //                ++cj;
            //                break;
            //            case TouchType.ConjoinedMove:
            //                ++cM;
            //                break;
            //            case TouchType.Pinch:
            //                ++p;
            //                break;
            //            case TouchType.Single:
            //                ++s;
            //            break;
            //            case TouchType.SingleMove:
            //                ++sM;
            //            break;
            //        }
            //    }

            //    Debug.LogFormat("Cj:{0}  Cm:{1}  P:{2} s:{3} sM:{4}", cj, cM, p, s, sM);
            //}

            #endregion
        }

        void TryRegisterTap(ExtendedTouch et)
        {
            foreach (var t in et.touches)
            {
                if (Time.time - prevTouchDatasP[t.fingerId].time > tapTimeout || Time.time - prevTouchDatasP[t.fingerId].time<0.025f)
                    continue;

                if (Vector2.Distance(t.position, prevTouchDatasP[t.fingerId].position)/Dpi > tapMaxRadius)
                    continue;

                if (Time.time - prevTouchDatasP[t.fingerId].lastTapTime < doubleTapTimeout && Time.time - prevTouchDatasP[t.fingerId].lastTapTime > 0.025f)
                {
                    if (EDoubleTap != null)
                        EDoubleTap(t.position);
                }
                else
                {
                    if (ETap != null)
                        ETap(t.position);
                }

                prevTouchDatasP[t.fingerId].lastTapTime = Time.time;
            }
        }

        private void RegisterPrevTouches(List<Touch> touches)
        {
            foreach (var t in touches)
            {
                if (t.phase != TouchPhase.Began)
                    continue;

                prevTouchDatasP[t.fingerId].position = t.position;
                prevTouchDatasP[t.fingerId].time = Time.time;
            }
        }

        private bool PassTouchProtection(int index, int id)
        {
            if (!touchProtectorsP[index].protect)
                return true;

            if (touchProtectorsP[index].hostId == id)
                return true;

            return false;
        }

        public void ProtectToutch(int extId)
        {
            ExtendedTouch exTouch;
            GetTouch(extId, out exTouch);

            foreach (var t in exTouch.touches)
            {
                touchProtectorsP[t.fingerId].protect = true;
                touchProtectorsP[t.fingerId].hostId = extId;
            }
        }

        private void RemoveProtection(ExtendedTouch extTouch)
        {
            foreach (var t in extTouch.touches)
            {
                touchProtectorsP[t.fingerId].protect = false;
            }
        }

        private void ApplyDeadZones(List<Touch> touches)
        {
            for (int i = 0; i < touches.Count; i++)
            {
                Touch t = touches[i];
                if (!(Math.Abs(t.deltaTime) > Mathf.Epsilon))
                    continue;

                if (t.phase != TouchPhase.Moved)
                    continue;

                var delta = t.deltaPosition/Dpi/t.deltaTime;
                if (delta.magnitude<deadZone)
                {
                    t.phase = TouchPhase.Stationary;
                    touches[i] = t;
                }
            }
        }

        private void MarkAsUsed(List<ExtendedTouch> extTouches, List<Touch> touches)
        {
            foreach (var ot in extTouches)
            {
                foreach (var t in ot.touches)
                {
                    for (int i = 0; i < touches.Count; i++)
                    {
                        var ct = touches[i];
                        if (ct.fingerId == t.fingerId)
                        {
                            ct.phase = TouchPhase.Canceled;
                            touches[i] = ct;
                        }
                    }
                }
            }
        }

        private void FindEnded(List<ExtendedTouch> newTouches)
        {
            for (int i = 0; i < extendedTouches.Count; i++)
            {
                var ot = extendedTouches[i];
                if (newTouches.All(t => !t.Equals(ot)))
                {
                    ot.phase = TouchPhase.Ended;
                    newTouches.Add(ot);
                    RemoveProtection(ot);
                    TryRegisterTap(ot);
                }
            }

            extendedTouches = newTouches;
        }

        private void FindSingle(List<Touch> touches, List<ExtendedTouch> resultTouches, List<ExtendedTouch> oldTouches)
        {
            for (int i = 0; i < touches.Count; i++)
            {
                var t = touches[i];
                if (t.phase == TouchPhase.Canceled || t.phase == TouchPhase.Ended)
                    continue;

                var hash = (t.fingerId + 1)*3 + TouchType.Single.GetHashCode();

                if (!PassTouchProtection(t.fingerId, hash))
                    continue;

                var ths = new List<Touch>() {t};
                var delta = t.deltaPosition;
                var index = oldTouches.FindIndex(ot => ot.GetHashCode() == hash);
                if (index == -1)
                {
                    if(!BlockedByGUI(t))
                        resultTouches.Add(new ExtendedTouch(TouchType.Single, TouchPhase.Began, ths, delta));
                }
                else
                {
                    var rt = oldTouches[index];
                    rt.phase = t.phase;
                    rt.Delta = t.deltaPosition;
                    rt.touches[0] = t;
                    resultTouches.Add(rt);
                }

                t.phase = TouchPhase.Canceled;
                touches[i] = t;
            }
        }

        private void FindPinch(List<Touch> touches, List<ExtendedTouch> resultTouches, List<ExtendedTouch> oldTouches)
        {
            for (int first = 0; first < touches.Count; first++)
            {
                var fst = touches[first];
                for (int second = 0; second < touches.Count; second++)
                {
                    var snd = touches[second];

                    if (!NeedToCheck(fst, snd))
                        continue;

                    var hash = (fst.fingerId + 1)*3 + (snd.fingerId + 1)*3 + TouchType.Pinch.GetHashCode();

                    if (!PassTouchProtection(fst.fingerId, hash) || !PassTouchProtection(snd.fingerId, hash))
                        continue;


                    float startDist = Vector3.Distance(fst.position, snd.position);
                    float newDist = Vector3.Distance(fst.position+fst.deltaPosition, snd.position+snd.deltaPosition);
                    float delta = startDist - newDist;
                    Vector3 deltaVector = Vector3.one * delta;

                   

                    Vector3 fdN = fst.deltaPosition.normalized;
                    Vector3 sdN = snd.deltaPosition.normalized;
                    float deltaDot = Vector2.Dot(fdN, sdN);

                    Vector2 center = snd.position + (snd.position - fst.position) / 2;

                    float first2center = 0;
                    float second2center = 0;

                    if (fst.deltaPosition.sqrMagnitude/Dpi > 0.01f)
                    {
                        first2center = Math.Abs(Vector3.Dot(fdN, (fst.position - center).normalized));
                    }

                    if (snd.deltaPosition.sqrMagnitude / Dpi > 0.01f)
                    {
                        second2center = Math.Abs(Vector3.Dot(sdN, (snd.position - center).normalized));
                    }

                    bool canCreateNew = Math.Abs(deltaDot) > 0 && Math.Max(first2center, second2center) > 0.9 && !(BlockedByGUI(fst) || BlockedByGUI(snd));
                    
                    Action<int, int> clearUsed = (fInd,sInd) =>
                    {
                        fst.phase = TouchPhase.Canceled;
                        snd.phase = TouchPhase.Canceled;

                        touches[fInd] = fst;
                        touches[sInd] = snd;
                    };
                    
                    var index = oldTouches.FindIndex(ot => ot.GetHashCode() == hash);
                    if (index != -1)
                    {
                        if (BlockedByGUI(fst) || BlockedByGUI(snd))
                            return;

                        TouchPhase phase = deltaDot < 0.5f && Math.Max(first2center, second2center) > 0.5
                            ? TouchPhase.Moved
                            : TouchPhase.Stationary;

                        if (phase == TouchPhase.Moved)
                        {
                            var rt = oldTouches[index];
                            rt.Delta = deltaVector;
                            rt.touches[0] = fst;
                            rt.touches[1] = snd;
                            rt.phase = TouchPhase.Moved;
                            resultTouches.Add(rt);
                            clearUsed(first, second);
                        }
                        else
                        {
                            var rt = oldTouches[index];
                            rt.Delta = deltaVector;
                            rt.touches[0] = fst;
                            rt.touches[1] = snd;
                            rt.phase = TouchPhase.Stationary;
                            resultTouches.Add(rt);
                            clearUsed(first, second);
                        }
                    }
                    else if(canCreateNew)
                    {
                        resultTouches.Add(new ExtendedTouch(TouchType.Pinch, TouchPhase.Began, new List<Touch>() { fst, snd }, deltaVector));
                        clearUsed(first, second);
                    }
                    
                }
            }
        }

        private void FindTwist(List<Touch> touches, List<ExtendedTouch> resultTouches, List<ExtendedTouch> oldTouches, bool disableAfterSuccess)
        {
            for (int first = 0; first < touches.Count; first++)
            {
                var fst = touches[first];
                for (int second = 0; second < touches.Count; second++)
                {
                    var snd = touches[second];

                    if (UsedInThisOldToutches(resultTouches, TouchType.Twist, fst, snd))
                        return;

                    if (!NeedToCheck(fst, snd))
                        continue;

                    var hash = (fst.fingerId + 1)*3 + (snd.fingerId + 1)*3 + TouchType.Twist.GetHashCode();

                    if (!PassTouchProtection(fst.fingerId, hash) || !PassTouchProtection(snd.fingerId, hash))
                        continue;


                    Vector2 fN = fst.deltaPosition.normalized;
                    Vector2 sN = snd.deltaPosition.normalized;
                    float deltaDot = Vector2.Dot(fN, sN);

                    Vector2 center = snd.position + (snd.position - fst.position) / 2;
                    
                    Vector2 f2c = (fst.position - center).normalized;
                    Vector2 s2c = (snd.position - center).normalized;
                    
                    Vector2 fTan = new Vector2(-f2c.y, f2c.x);
                    Vector2 sTan = new Vector2(-s2c.y, f2c.x);

                    bool isTangent = Mathf.Abs(Vector2.Dot(fN, fTan))>0.95f || Mathf.Abs(Vector2.Dot(sN, sTan)) > 0.95f;

                    float delta = 0;
                    if (fst.deltaPosition.sqrMagnitude < 0.01f)
                    {
                        delta = TwistDelta(snd, fst.position);
                    }
                    else if (snd.deltaPosition.sqrMagnitude < 0.01f)
                    {
                        delta = TwistDelta(fst, snd.position);
                    }
                    else
                    {
                        delta = TwistDelta(fst, center)*2;
                    }

                    bool canCreateNew = deltaDot < 0 && isTangent && !(BlockedByGUI(fst) || BlockedByGUI(snd));

                    Action<int, int> clearUsed = (fInd, sInd) =>{ };
                    if (disableAfterSuccess)
                    {
                        clearUsed = (fInd, sInd) =>
                        {
                            fst.phase = TouchPhase.Canceled;
                            snd.phase = TouchPhase.Canceled;

                            touches[fInd] = fst;
                            touches[sInd] = snd;
                        };
                    }

                    var index = oldTouches.FindIndex(ot => ot.GetHashCode() == hash);
                    if (index != -1)
                    {
                        if (BlockedByGUI(fst) || BlockedByGUI(snd))
                            return;

                        TouchPhase phase = Mathf.Abs(delta) > 0.001f
                            ? TouchPhase.Moved
                            : TouchPhase.Stationary;

                        if (phase == TouchPhase.Moved)
                        {
                            var rt = oldTouches[index];
                            rt.Delta = new Vector2(delta, delta);/////////////////////////
                            rt.touches[0] = fst;
                            rt.touches[1] = snd;
                            rt.phase = TouchPhase.Moved;
                            resultTouches.Add(rt);
                            clearUsed(first, second);
                        }
                        else
                        {
                            var rt = oldTouches[index];
                            rt.Delta = Vector2.zero;/////////////////////////
                            rt.touches[0] = fst;
                            rt.touches[1] = snd;
                            rt.phase = TouchPhase.Stationary;
                            resultTouches.Add(rt);
                            clearUsed(first, second);
                        }
                    }
                    else if (canCreateNew)
                    {
                        resultTouches.Add(new ExtendedTouch(TouchType.Twist, TouchPhase.Began, new List<Touch>() { fst, snd }, Vector2.one));
                        clearUsed(first, second);
                    }

                }
            }
        }

        private void FindCojoined(List<Touch> touches, List<ExtendedTouch> resultTouches, List<ExtendedTouch> oldTouches)
        {
            for (int first = 0; first < touches.Count; first++)
            {
                var fst = touches[first];
                for (int second = 0; second < touches.Count; second++)
                {
                    var snd = touches[second];

                    if (!NeedToCheck(fst, snd))
                        continue;

                    var hash = (fst.fingerId + 1)*3 + (snd.fingerId + 1)*3 + TouchType.Conjoined.GetHashCode();

                    if (!PassTouchProtection(fst.fingerId, hash) || !PassTouchProtection(snd.fingerId, hash))
                        continue;

                    if (UsedInOtherOldToutches(oldTouches, TouchType.Conjoined, fst, snd))
                        return;

                    if (Vector2.Distance(fst.position, snd.position) / Dpi < cojointDistance/25.4f)
                    {
                        //8text.text = string.Format("dpi{2}\ndInI:{0}\nraw:{1}", Vector2.Distance(fst.position, snd.position)/ Dpi, Vector2.Distance(touch.position, snd.position), Screen.dpi);

                        var ths = new List<Touch>() { fst, snd };
                        var delta = (fst.deltaPosition + snd.deltaPosition) / 2f;


                        var index = oldTouches.FindIndex(ot => ot.GetHashCode() == hash);
                        if (index == -1)
                        {
                            if(BlockedByGUI(fst)||BlockedByGUI(snd))
                                return;
                            resultTouches.Add(new ExtendedTouch(TouchType.Conjoined, TouchPhase.Began, ths, delta));
                        }
                        else
                        {
                            var rt = oldTouches[index];

                            rt.phase = (fst.phase == TouchPhase.Moved || snd.phase == TouchPhase.Moved)
                                ? TouchPhase.Moved
                                : TouchPhase.Stationary;

                            rt.Delta = fst.deltaPosition.sqrMagnitude > snd.deltaPosition.sqrMagnitude
                                ? fst.deltaPosition
                                : snd.deltaPosition;
                            rt.touches[0] = fst;
                            rt.touches[1] = snd;

                            resultTouches.Add(rt);
                        }

                        fst.phase = TouchPhase.Canceled;
                        snd.phase = TouchPhase.Canceled;

                        touches[first] = fst;
                        touches[second] = snd;
                    }
                }
            }
        }


        private float TwistDelta(Touch touch, Vector2 center)
        {
            Vector2 oldPos = touch.position - touch.deltaPosition;
            Vector2 newVec = touch.position - center;
            Vector2 oldVec = oldPos - center;
            float result = Mathf.DeltaAngle(Mathf.Atan2(newVec.x, newVec.y) * Mathf.Rad2Deg, Mathf.Atan2(oldVec.x, oldVec.y) * Mathf.Rad2Deg);
            return result;
        }
        
        bool NeedToCheck(Touch fst, Touch snd)
        {
            if (fst.fingerId == snd.fingerId)
                return false;
            if (fst.phase == TouchPhase.Canceled || fst.phase == TouchPhase.Ended || snd.phase == TouchPhase.Canceled || snd.phase == TouchPhase.Ended)
                return false;

            //if (resultTouches.Any(rt => rt.touches.Any(t => t.fingerId == snd.fingerId)))
            //    return false;

            return true;
        }

        bool UsedInOtherOldToutches(List<ExtendedTouch> oldTouches, TouchType type, params Touch[] touches)
        {
            foreach (var ot in oldTouches)
            {
                if (ot.type == type || ot.type == TouchType.Single)
                    continue;

                foreach (var t in ot.touches)
                {
                    if (touches.Any(ct => ct.fingerId == t.fingerId))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        bool UsedInThisOldToutches(List<ExtendedTouch> oldTouches, TouchType type, params Touch[] touches)
        {
            foreach (var ot in oldTouches)
            {
                if (ot.type != type)
                    continue;

                foreach (var t in ot.touches)
                {
                    if (touches.Any(ct => ct.fingerId == t.fingerId))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        public bool GetTouch(int id, out ExtendedTouch result)
        {
            int index = extendedTouches.FindIndex(t => t.GetHashCode() == id);
            if (index != -1)
            {
                result = extendedTouches[index];
                return true;
            }

            result = new ExtendedTouch();
            return false;
        }

        public bool GetFirstTouchWithPhase(TouchPhase phase, out ExtendedTouch result)
        {
            int index = extendedTouches.FindIndex(t => t.phase == phase);
            if (index!=-1)
            {
                result = extendedTouches[index];
                return true;
            }

            result = new ExtendedTouch();
            return false;
        }

        public bool GetFirstTouchWithType(TouchType type, out ExtendedTouch result)
        {
            int index = extendedTouches.FindIndex(t => t.type == type);
            if (index != -1)
            {
                result = extendedTouches[index];
                return true;
            }

            result = new ExtendedTouch();
            return false;
        }


        public List<ExtendedTouch> GetTouches(TouchPhase phase)
        {
            return extendedTouches.Where(t => t.phase == phase).ToList();
        }

        public List<ExtendedTouch> GetTouches(TouchType type)
        {
            return extendedTouches.Where(t => t.type == type).ToList();
        }

        public void OnGUI()
        {
            foreach (var et in extendedTouches)
            {
                switch (et.type)
                {
                    case TouchType.Conjoined:
                        DrawTouch(et.Center, 2.5f, new Color(0.25f, 1, 0.25f, 0.25f));
                        break;

                    case TouchType.Single:
                        DrawTouch(et.Center, 1.4f, new Color(0.25f, 0.5f, 1, 0.25f));
                        break;

                    case TouchType.Pinch:
                    case TouchType.Twist:
                        DrawTouch(et.touches[0].position, 1.4f, new Color(0.25f, 0.25f, 1, 0.25f));
                        DrawTouch(et.touches[1].position, 1.4f, new Color(1, 0.25f, 0.25f, 0.25f));
                        DrawTouch(et.Center, 0.2f, new Color(0.25f, 0.25f, 0.25f, 0.5f));
                        break;
                }
            }
        }

        void DrawTouch(Vector2 pos, float sizeMtl, Color color)
        {
            if(!touchCircle)
                return;

            var tPos = pos;
            tPos.y = Screen.height - tPos.y;
            var sizeX = Dpi/2*Camera.main.aspect*sizeMtl;
            var sizeY = Dpi / 2 * sizeMtl;
            sizeX = (int)(sizeX / Camera.main.aspect);
            Color oldColor = GUI.color;
            GUI.color = color;
            GUI.DrawTexture(new Rect(tPos.x - sizeX / 2, tPos.y - sizeY / 2, sizeX, sizeY), touchCircle);
            GUI.color = oldColor;
        }

        public bool BlockedByGUI(Touch touch)
        {

            PointerEventData eventDataCurrentPosition = new PointerEventData(EventSystem.current);
            eventDataCurrentPosition.position = touch.position;
            List<RaycastResult> results = new List<RaycastResult>();
            EventSystem.current.RaycastAll(eventDataCurrentPosition, results);
            return results.Any(r => (1 << r.gameObject.layer | uiMask.value) == uiMask);
        }
    }

    [System.Serializable]
    public struct ExtendedTouch
    {
        public ExtendedTouch(TouchType type, TouchPhase phase,List<Touch> touches, Vector2 delta)
        {
            this.type = type;
            this.touches = touches;
            this.Delta = delta;
            this.phase = phase;
            id = -1;
            id = GetHashCode();
        }

        public TouchType type;
        public TouchPhase phase;
        public int id;

        public List<Touch> touches;
        public Vector2 Delta;

        public override int GetHashCode()
        {
            {
                switch (type)
                {
                    case TouchType.Conjoined:
                    case TouchType.Pinch:
                    case TouchType.Twist:
                        return type.GetHashCode() + (touches[0].fingerId+1)*3 + (touches[1].fingerId+1)*3;
                    case TouchType.Single:
                        return type.GetHashCode() + (touches[0].fingerId+1)*3;
                }
                return 0;
            }
        }

        public override bool Equals(object other)
        {
            if (other == null)
                return false;

            return GetHashCode() == ((ExtendedTouch)other).GetHashCode();
        }

        public Vector2 Center
        {
            get
            {
                Vector2 c = new Vector2();

                if (touches.Count == 0)
                    return c;

                foreach (var t in touches)
                    c += t.position;

                return c/touches.Count;
            }
        }

        public float DeltaTime
        {
            get
            {
                float summ = 0;
                for (int i = 0; i < touches.Count; i++)
                    summ += touches[i].deltaTime;

                return summ/touches.Count;
            }
        }
    }

    public enum TouchType
    {
        Conjoined = 100,
        Pinch = 200,
        Single = 300,
        Twist = 400,
    }

    [System.Serializable]
    public struct TouchProtector
    {
        public bool protect;
        public int hostId;
    }

    [System.Serializable]
    public struct PrevTouchData
    {
        public Vector2 position;
        public float time;
        public float lastTapTime;
    }
}

namespace UnityEngine.EventSystems
{
}