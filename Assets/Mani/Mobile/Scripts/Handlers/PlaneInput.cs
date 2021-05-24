using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public class PlaneInput : DragInput
    {
        private Vector3 beginPointP;
        private Vector3 beginShiftP;

        public override void OnDrag(PointerEventData eventData)
        {
            if (!handle)
                return;

            Plane plane = new Plane(transform.right, beginPointP);
            float dist;
            Ray ray = Camera.main.ScreenPointToRay(eventData.position);
            if (!plane.Raycast(ray, out dist))
                return;

            Vector3 newPos = ray.GetPoint(dist);
            handle.OnDrag(newPos-beginShiftP, true);
        }

        public override void OnBeginDrag(PointerEventData eventData)
        {
            beginPointP = transform.parent.position;

            Plane plane = new Plane(transform.right, beginPointP);
            float dist;
            Ray ray = Camera.main.ScreenPointToRay(eventData.position);
            if (!plane.Raycast(ray, out dist))
                return;

            beginShiftP = ray.GetPoint(dist) - beginPointP;

            base.OnBeginDrag(eventData);
        }
    }
}