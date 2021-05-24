using System;
using System.Collections.Generic;
using System.Linq;
using ManipulationInstruments.Handles;
using UnityEngine;

namespace ManipulationInstruments
{
    public class ManipulatorsController : MonoBehaviour
    {
        
        public event Action EOnBeginDrag;

        public static ManipulatorsController Instance
        {
            get { return instanceS; }
        }
        private static ManipulatorsController instanceS;

        void Awake()
        {
            if (!instanceS)
                instanceS = this;
            else
                Debug.LogError("Must be only one");
        }

		public GameObject CoilSettings;
        public GameObject AxesRotation;

        private HandleBase currentInstrument;
        public HandleBase[] manipulatorsPrefabs;
        private List<HandleBase> manipulatorsP = new List<HandleBase>();
        private List<GameObject> lastSelectionP;

        public Coil prevCoil;

        private void OnUpdateSelection(List<GameObject> secection)
        {
            lastSelectionP = secection;
            currentInstrument.UpdateSelection(secection);

            if (secection.Count == 0)
                CoilSettings.SetActive(false);
            else
            {
                if (prevCoil != null) 
                    SetSettings();
                
                ReadSettings(secection);
                CoilSettings.SetActive(true);
            }
        }

        private void ReadSettings(List<GameObject> selection)
        {
            var Coil = selection[0].GetComponent<Coil>();
            prevCoil = Coil;

            Coil.SetUIValues();
        }

		private void SetSettings()
        {
            prevCoil.GetUIValues();
        }
        
        void Start()
        {
            MouseSelector.EOnSelectionCheange += OnUpdateSelection;

            foreach (var prefab in manipulatorsPrefabs)
            {
                HandleBase inst = Instantiate(prefab);
                inst.Deactivate();
                manipulatorsP.Add(inst);
            }

            currentInstrument = manipulatorsP.First(m => m.handleType == HandleType.Translate);
            currentInstrument.Activate();
        }

        void OnDestroy()
        {
            MouseSelector.EOnSelectionCheange -= OnUpdateSelection;
        }

        public static void UpdateSelection()
        {
            Instance.currentInstrument.UpdateSelection(Instance.lastSelectionP);
        }

        public static void CheangeInstrument(HandleType type)
        {
            Instance.currentInstrument.Deactivate();
            Instance.currentInstrument = instanceS.manipulatorsP.First(m => m.handleType == type);
            Instance.currentInstrument.UpdateSelection(Instance.lastSelectionP);
            Instance.currentInstrument.Activate();
        }

        public void OnBeginDrag()
        {
            if (EOnBeginDrag != null)
                EOnBeginDrag();
        }
    }
}