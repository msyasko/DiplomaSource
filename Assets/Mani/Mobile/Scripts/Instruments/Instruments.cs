using Gestures;
using UnityEngine;
using TouchType = Gestures.TouchType;

namespace Instruments
{
    public abstract class Instruments : MonoBehaviour
    {
        public LayerMask selectableMask;

        protected Mode instrumentModeP;

        protected GameObject GetSelectebleFromSingleTouch()
        {
            ExtendedTouch singleTouch;
            if (!TouchInputManager.Instance.GetFirstTouchWithType(TouchType.Single, out singleTouch))
                return null;
            if (singleTouch.phase != TouchPhase.Began)
                return null;

            RaycastHit hit;
            if (Physics.Raycast(Camera.main.ScreenPointToRay(singleTouch.Center), out hit, Mathf.Infinity, selectableMask))
            {
                Debug.DrawRay(hit.point, Vector3.up * 3, Color.magenta);
                return hit.transform.root.gameObject;
            }

            return null;
        }
        
        protected GameObject GetSelectebleFromCojoinedTouch()
        {
            ExtendedTouch cojTouch;
            if (!TouchInputManager.Instance.GetFirstTouchWithType(TouchType.Conjoined, out cojTouch))
                return null;
            if (cojTouch.phase != TouchPhase.Began)
                return null;

            RaycastHit hit;
            if (Physics.Raycast(Camera.main.ScreenPointToRay(cojTouch.Center), out hit, Mathf.Infinity, selectableMask))
            {
                Debug.DrawRay(hit.point, Vector3.up * 3, Color.magenta);
                return hit.transform.root.gameObject;
            }

            return null;
        }

        protected virtual void Awake()
        {
            InstrumentsHub.EModeDisablad += InstrumentDisablad;
            InstrumentsHub.EModeEnabled += InstrumentEnabled;
            InitType();
        }

        void OnDestroy()
        {
            InstrumentsHub.EModeDisablad -= InstrumentDisablad;
            InstrumentsHub.EModeEnabled -= InstrumentEnabled;
        }

        void InstrumentDisablad(Mode mode)
        {
            if (instrumentModeP == mode)
                enabled = false;
        }

        void InstrumentEnabled(Mode mode)
        {
            if (instrumentModeP == mode)
                enabled = true;
        }

        protected abstract void InitType();
    }
}