using System.Collections.Generic;
using UnityEngine;

namespace ManipulationInstruments.Handles
{
    public class RotationHandle : HandleBase
    {
        public List<Quaternion> stardLocalRotations = new List<Quaternion>();

        public override void OnDrag(Vector3 vector, bool dorectMode)
        {
            transform.RotateAround(transform.position, vector, vector.magnitude);
            if (!IsMultiselect)
            {
                foreach (var o in selection)
                    o.transform.rotation = transform.rotation;
            }
            else
            {
                for (int i = 0; i < selection.Count; i++)
                {
                    selection[i].transform.position = transform.position + transform.TransformDirection(stardLocalPoses[i]);
                    selection[i].transform.rotation = transform.rotation * stardLocalRotations[i];
                }
            }
        }

        public override void UpdateSelection(List<GameObject> selection)
        {
            if(selection == null)
                return;
            base.UpdateSelection(selection);

            stardLocalRotations.Clear();

            foreach (var o in selection)
            {
                stardLocalRotations.Add(Quaternion.Inverse(transform.rotation)*o.transform.rotation);
            }
        }

        //public override void OnBeginDrag()
        //{
        //    base.OnBeginDrag();
        //    Debug.Log("RotateBegin");
        //}

        //public override void OnEndDrag()
        //{
        //    base.OnEndDrag();
        //    Debug.Log("RotateEnd");
        //}
    }
}