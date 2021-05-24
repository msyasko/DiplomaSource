using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public class RingInput : DragInput
    {
        public override void OnDrag(PointerEventData eventData)
        {
            if (!handle)
                return;

            Vector3 delta = eventData.delta;
            delta /= MinScreenSize;
            Vector3 localDelta = transform.InverseTransformDirection(Camera.main.transform.TransformDirection(delta));
            Vector3 cameraForward = Camera.main.transform.forward;
            Vector3 camUpCross = Vector3.Cross(cameraForward, transform.forward);
            float direction = Vector3.Dot(transform.TransformDirection(localDelta).normalized, camUpCross);
            direction = Mathf.Sign(direction);
            handle.OnDrag(transform.forward * localDelta.magnitude * 90 * direction, false);
        }
    }
}