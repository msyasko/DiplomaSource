using System.Collections.Generic;
using UnityEngine;

namespace ManipulationInstruments.Handles
{
    public abstract class HandleBase : MonoBehaviour
    {
        private SceneStateHistory undoRedoControllerP;

        private HistoryState_TransformMultiselect snapshot;

        public HandleType handleType;
        public bool useWorldSpace;
        public float gizmoSizeMtl = 5;

        protected List<GameObject> selection = new List<GameObject>();
        protected List<Vector3> stardLocalPoses = new List<Vector3>();

        public static bool InDragMode { get; private set; }
        
        public abstract void OnDrag(Vector3 vector, bool dorectMode);

        public bool IsMultiselect { get { return selection.Count > 1; } }

        void Start()
        {
            undoRedoControllerP = UnitySingleton<SceneStateHistory>.Instance;
        }

        public virtual void OnBeginDrag()
        {
            InDragMode = true;
            ManipulatorsController.Instance.OnBeginDrag();
            SnapOld();
        }

        private void SnapOld()
        {
            snapshot = new HistoryState_TransformMultiselect();
            
            foreach (var o in selection)
            {
                snapshot.oldSnap.Add(new HistoryState_TransformMultiselect.TransformSnapshot()
                {
                    transform = o.transform,
                    position = o.transform.position,
                    rotation = o.transform.rotation,
                    scale = o.transform.localScale
                });
            }
        }

        public virtual void OnEndDrag()
        {
            InDragMode = false;

            SnapNew();
            undoRedoControllerP.AddState(snapshot);
        }

        private void SnapNew()
        {
            foreach (var o in selection)
            {
                snapshot.newSnap.Add(new HistoryState_TransformMultiselect.TransformSnapshot()
                {
                    transform = o.transform,
                    position = o.transform.position,
                    rotation = o.transform.rotation,
                    scale = o.transform.localScale
                });
            }
        }

        public void SwithSpace()
        {
            useWorldSpace = !useWorldSpace;
            if (selection.Count == 0)
                return;
            transform.rotation = useWorldSpace ? Quaternion.identity : selection[0].transform.rotation;
        }

        public virtual void UpdateSelection(List<GameObject> selection)
        {
            if(selection == null)
                return;
            this.selection = selection;
            if (selection.Count == 0)
                Deactivate();
            else
            {
                gameObject.SetActive(true);
                PlaceInSelectionCenter();
                AliginWithSelection(selection);
                UpdateStartLocalPoses(selection);
            }
        }

        private void UpdateStartLocalPoses(List<GameObject> selection)
        {
            stardLocalPoses.Clear();
            foreach (var o in selection)
                stardLocalPoses.Add(transform.InverseTransformDirection(o.transform.position - transform.position));
        }

        private void AliginWithSelection(List<GameObject> selection)
        {
            if (selection.Count == 1)
                transform.rotation = useWorldSpace ? Quaternion.identity : selection[0].transform.rotation;
            else
                transform.rotation = Quaternion.identity;
        }

        public void Deactivate()
        {
            gameObject.SetActive(false);
        }

        public void Activate()
        {
            if (selection.Count == 0)
                return;

            gameObject.SetActive(true);
            PlaceInSelectionCenter();
            AliginWithSelection(selection);
            UpdateStartLocalPoses(selection);
        }

        private void PlaceInSelectionCenter()
        {
            Vector3 center = Vector3.zero;
            foreach (var o in selection)
                center += o.transform.position;
            center /= selection.Count;
            transform.position = center;
        }

        void LateUpdate()
        {
            float size = GetGizmoSize(transform.position);
            transform.localScale = Vector3.one * size * gizmoSizeMtl;
        }

        public static float GetGizmoSize(Vector3 position)
        {
            Camera current = Camera.main;

            if (current)
            {
                Transform transform = current.transform;
                Vector3 position2 = transform.position;
                float z = Vector3.Dot(position - position2, transform.TransformDirection(new Vector3(0f, 0f, 1f)));
                Vector3 a = current.WorldToScreenPoint(position2 + transform.TransformDirection(new Vector3(0f, 0f, z)));
                Vector3 b = current.WorldToScreenPoint(position2 + transform.TransformDirection(new Vector3(1f, 0f, z)));
                float magnitude = (a - b).magnitude;
                return 80f / Mathf.Max(magnitude, 0.0001f);
            }

            return 20f;
        }

        //void Update()
        //{
        //    for (int i = 0; i < selection.Count; i++)
        //    {
        //        Debug.DrawLine(transform.position, transform.position +
        //                                           transform.TransformDirection(stardLocalPosesP[i]), Color.yellow);

        //        Debug.DrawRay(selection[i].transform.position, Vector3.up, Color.blue);
        //    }
        //}
    }

    public enum HandleType
    {
        Translate, Scale, Rotate
    }
}