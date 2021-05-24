using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public class UniformInput : DragInput
    {
        private Vector3 directionToPointer;

        public override void OnDrag(PointerEventData eventData)
        {
            if (!handle)
                return;

            Vector3 delta = eventData.delta;
            delta /= MinScreenSize;

            
            float direction = Vector3.Dot(directionToPointer, delta);
            direction = Mathf.Sign(direction);
            handle.OnDrag(Vector3.one*direction*delta.magnitude, true);
        }

        public override void OnBeginDrag(PointerEventData eventData)
        {
            Vector3 positionRelCamera = Camera.main.WorldToScreenPoint(transform.position);
            positionRelCamera.z = 0;
            Vector3 pointerPos = eventData.position;
            directionToPointer = (pointerPos - positionRelCamera).normalized;
            base.OnBeginDrag(eventData);
        }
    }
}