using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;

namespace Mobile.UI
{
    public class UGUIButtonEventsHandler : MonoBehaviour, IPointerDownHandler, IPointerUpHandler, IPointerClickHandler
    {
        public UnityEvent onDownEvent;
        public UnityEvent onUpEvent;
        public UnityEvent onClickEvent;
        public UnityEvent onDoubleClickEvent;

        private float lastClickTimeP;
        private float doubleClickTimeoutP = 0.3f; 

        public void OnPointerDown(PointerEventData eventData)
        {
            onDownEvent.Invoke();
        }

        public void OnPointerUp(PointerEventData eventData)
        {
            onUpEvent.Invoke();
        }

        public void OnPointerClick(PointerEventData eventData)
        {
            if (Time.time - lastClickTimeP < doubleClickTimeoutP)
                OnDoubleClick();
            else
            {
                lastClickTimeP = Time.time;
                onClickEvent.Invoke();
            }
        }

        private void OnDoubleClick()
        {
            onDoubleClickEvent.Invoke();
        }
    }
}