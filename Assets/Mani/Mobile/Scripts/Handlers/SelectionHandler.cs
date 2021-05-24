using System.Collections.Generic;
using UnityEngine;

namespace ManipulationInstruments
{
    public class SelectionHandler:MonoBehaviour
    {
        public GameObject manipulators;
        

        void Start()
        {
            MouseSelector.EOnSelectionCheange += OnSelect;

            HideAll();
        }

        void OnDestroy()
        {
            MouseSelector.EOnSelectionCheange -= OnSelect;
        }

        private void HideAll()
        {
            manipulators.SetActive(false);
        }

        private void OnSelect(List<GameObject> selection)
        {
            manipulators.SetActive(selection != null && selection.Count > 0);
        }
    }
}