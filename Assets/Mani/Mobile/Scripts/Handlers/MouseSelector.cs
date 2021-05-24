using System;
using System.Collections.Generic;
using System.Linq;
using Gestures;
using Instruments;
using UnityEngine;
using UnityEngine.EventSystems;
using ManipulationInstruments.Handles;
using TouchType = Gestures.TouchType;

namespace ManipulationInstruments
{
    public class MouseSelector : MonoBehaviour
    {
        public static MouseSelector Instance { get { return instanceS; } }
        private static MouseSelector instanceS;
        
        public bool multiselectMode;

        public static event Action<List<GameObject>> EOnSelectionCheange; 

        public LayerMask uiMask;
        public LayerMask selectableLayers;

        public List<GameObject> selection = new List<GameObject>();
        private SceneStateHistory undoRedoControllerP;

        void Awake()
        {
            if (!instanceS)
                instanceS = this;

            StateController.ESetState += OnSetState;
        }

        void OnDestroy()
        {
            StateController.ESetState -= OnSetState;
        }

        private void OnSetState(SceneState state)
        {
            if (state != SceneState.Setup)
                SelectSingle(null);
        }

        void Start()
        {
            undoRedoControllerP = UnitySingleton<SceneStateHistory>.Instance;
        }

        public void SelectSingle(GameObject go)
        {
            SelectCommand selectCommand = new SelectCommand(){oldOjObjects = new List<GameObject>(selection)};

            if (!multiselectMode)
            {
                selection.Clear();
                if (go != null)
                    selection.Add(go);
                OnSecectionSet();
                selectCommand.newOjObjects = new List<GameObject>(selection);
                undoRedoControllerP.AddState(selectCommand);
            }
            else
            {
                if (go!=null&&!selection.Contains(go))
                {
                    selection.Add(go);
                    OnSecectionSet();
                    selectCommand.newOjObjects = new List<GameObject>(selection);
                    undoRedoControllerP.AddState(selectCommand);
                }
            }
        }

        private void OnSecectionSet()
        {
            if (EOnSelectionCheange != null)
                EOnSelectionCheange(selection);
        }

        void Update()
        {
            if (HandleBase.InDragMode || InstrumentsHub.Instance.Mode != Mode.None)
                return;

            ExtendedTouch single;
            if (TouchInputManager.Instance.GetFirstTouchWithType(TouchType.Single, out single)
                &&
                single.phase==TouchPhase.Began)
            {
                Raycast(single.Center);
            }
            else if (Input.GetMouseButtonDown(0))
            {
                Raycast(new Vector2(Input.mousePosition.x, Input.mousePosition.y));
            }
        }

        void Raycast(Vector2 pos)
        {
            if (StateController.State != SceneState.Setup)
                return;

            Ray ray = Camera.main.ScreenPointToRay(pos);
            Debug.DrawRay(ray.origin, ray.direction * 100, Color.green, 1);

            if (IsPointerOverUIObject(pos))
                return;

            RaycastHit hit;
            if (Physics.Raycast(ray, out hit, Mathf.Infinity, selectableLayers))
                SelectSingle(GetRelativeRoot(hit.transform));
            else
                SelectSingle(null);
        }

        GameObject GetRelativeRoot(Transform tr)
        {
            GameObject result = tr.gameObject;

            Transform curent = tr;
            Transform root = tr.root;

            while (curent != root)
            {
                curent = curent.parent;
            }

            return result;
        }

        private bool IsPointerOverUIObject(Vector2 pos)
        {
            if (uiMask == -1)
                return true;

            PointerEventData eventDataCurrentPosition = new PointerEventData(EventSystem.current);
            eventDataCurrentPosition.position = pos;
            List<RaycastResult> results = new List<RaycastResult>();
            EventSystem.current.RaycastAll(eventDataCurrentPosition, results);
            var resBool = results.Any(r => (1 << r.gameObject.layer | uiMask.value) == uiMask);
            return resBool;
        }

        public void SetSingleMode()
        {
            multiselectMode = false;
        }

        public void SetMultipleMode()
        {
            multiselectMode = true;
        }

        public static void ForceSelect(List<GameObject> objects)
        {
            Instance.selection = objects;
            Instance.OnSecectionSet();
        }
    }

    public class SelectCommand : IHistoryState
    {
        public List<GameObject> oldOjObjects;
        public List<GameObject> newOjObjects;

        public void Dispose(){}

        public bool Execute()
        {
            MouseSelector.ForceSelect(newOjObjects);
            return true;
        }

        public bool Undo()
        {
            MouseSelector.ForceSelect(oldOjObjects);
            return true;
        }
    }
}