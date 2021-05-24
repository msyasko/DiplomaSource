using Instruments;
using UnityEngine;
using UnityEngine.Events;
using UnityEngine.EventSystems;

namespace GUIns
{
    public class UGUIButtonEventsHandlerForSidePanel : MonoBehaviour, IPointerDownHandler, IPointerUpHandler, IPointerClickHandler
    {
        public ModeEvent onDownEvent;
        public ModeEvent onUpEvent;
        public ModeEvent onClickEvent;
        public ModeEvent onDoubleClickEvent;

        private float lastClickTimeP;
        private float doubleClickTimeoutP = 0.3f;

        public Mode mode;

        public void OnPointerDown(PointerEventData eventData)
        {
            onDownEvent.Invoke(mode);
        }

        public void OnPointerUp(PointerEventData eventData)
        {
            onUpEvent.Invoke(mode);
        }

        public void OnPointerClick(PointerEventData eventData)
        {
            if (Time.time - lastClickTimeP < doubleClickTimeoutP)
                OnDoubleClick();
            else
            {
                lastClickTimeP = Time.time;
                onClickEvent.Invoke(mode);
            }
        }

        private void OnDoubleClick()
        {
            onDoubleClickEvent.Invoke(mode);
        }


        [System.Serializable]
        public class ModeEvent: UnityEvent<Mode>
        {
            
        }
    }
}