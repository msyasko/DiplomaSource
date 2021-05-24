using System.Collections.Generic;
using UnityEngine;

namespace ManipulationInstruments.Handles
{
    public class TranslateHandle : HandleBase
    {
        public override void OnDrag(Vector3 vector, bool dorectMode)
        {
            float distTocamera = Vector3.Distance(transform.position, Camera.main.transform.position);
            if (dorectMode)
                transform.position = vector;
            else
                transform.Translate(vector*distTocamera, Space.World);

            for (int i = 0; i < selection.Count; i++)
                selection[i].transform.position = transform.position+transform.TransformDirection(stardLocalPoses[i]);
        }

        //public override void OnBeginDrag()
        //{
        //    base.OnBeginDrag();
        //    Debug.Log("TranslateBegin");
        //}

        //public override void OnEndDrag()
        //{
        //    base.OnEndDrag();
        //    Debug.Log("TranslateEnd");
        //}
    }
}