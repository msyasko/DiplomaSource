using System;
using UnityEngine;

namespace Instruments
{
    public class InstrumentsHub : MonoBehaviour
    {
        public static InstrumentsHub Instance { get { return instanceP; } }
        private static InstrumentsHub instanceP;
        
        public static event Action<Mode> EModeDisablad;
        public static event Action<Mode> EModeEnabled;

        public InstrumentsHub()
        {
            Mode = Mode.None;
        }

        public Mode Mode { get; private set; }

        public void Activate(Mode newMode)
        {
            if (Mode != Mode.None)
                return;

            Mode = newMode;

            if (EModeEnabled != null)
                EModeEnabled(newMode);
        }

        public void Deactivate(Mode deactivatedMode)
        {
            if (Mode != deactivatedMode)
                return;

            Mode = Mode.None;

            if (EModeDisablad != null)
                EModeDisablad(deactivatedMode);
        }

        void Awake()
        {
            if (!Instance)
                instanceP = this;
        }
    }

    public enum Mode
    {
        None, Delete, Select, Freeze, Duplicate, BringToScene, PlaceOnTap
    }
}