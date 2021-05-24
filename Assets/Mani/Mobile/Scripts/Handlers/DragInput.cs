using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.Handles
{
    public abstract class DragInput : MonoBehaviour, IDragHandler, IBeginDragHandler, IEndDragHandler
    {
        public HandleBase handle;

        public abstract void OnDrag(PointerEventData eventData);
        
        public virtual void OnBeginDrag(PointerEventData eventData)
        {
            if (handle)
                handle.OnBeginDrag();
        }

        public virtual void OnEndDrag(PointerEventData eventData)
        {
            if (handle)
                handle.OnEndDrag();
        }

        protected int MinScreenSize { get { return Mathf.Min(Screen.width, Screen.height); } }
    }
}