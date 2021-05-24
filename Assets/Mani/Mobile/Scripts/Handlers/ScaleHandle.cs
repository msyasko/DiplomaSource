using UnityEngine;

namespace ManipulationInstruments.Handles
{
    public class ScaleHandle : HandleBase
    {
        public override void OnDrag(Vector3 vector, bool alternateMode)
        {
            for (int i = 0; i < selection.Count; i++)
            {
                if (!alternateMode)
                {
                    float distToCamera = Vector3.Distance(transform.position, Camera.main.transform.position);
                    selection[i].transform.localScale += transform.InverseTransformDirection(vector*distToCamera);
                }
                else
                    selection[i].transform.localScale += Vector3.Scale(selection[i].transform.localScale, vector);
            }
        }

        public override void OnBeginDrag()
        {
            base.OnBeginDrag();
            Debug.Log("ScaleBegin");
        }

        public override void OnEndDrag()
        {
            base.OnEndDrag();
            Debug.Log("ScaleEnd");
        }   
    }
}