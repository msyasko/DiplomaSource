using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public class ArrowInput : DragInput
    {
        public override void OnDrag(PointerEventData eventData)
        {
            if (!handle)
                return;

            Vector3 delta = eventData.delta;
            delta /= MinScreenSize;
            Vector3 localDelta = transform.InverseTransformDirection(Camera.main.transform.TransformDirection(delta));

            Vector3 direction = transform.TransformDirection(new Vector3(0, 0, localDelta.z)).normalized;
            handle.OnDrag(direction*localDelta.magnitude, false);
        }
    }
}
