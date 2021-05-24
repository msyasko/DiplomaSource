using ManipulationInstruments.Handles;
using UnityEngine;
using UnityEngine.EventSystems;

namespace ManipulationInstruments.UI
{
    public class SelectTrasnfomInstrumentButton:MonoBehaviour, IPointerClickHandler
    {
        public HandleType handleType;
        public GameObject deactevatedOnClick;

        public void OnPointerClick(PointerEventData eventData)
        {
            ManipulatorsController.CheangeInstrument(handleType);

            if (deactevatedOnClick != null)
                deactevatedOnClick.SetActive(false);
        }
    }
}