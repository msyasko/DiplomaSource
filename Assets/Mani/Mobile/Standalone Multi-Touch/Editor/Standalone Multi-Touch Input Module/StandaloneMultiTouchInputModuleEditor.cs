/* Written by Kaz Crowe */
/* StandaloneMultiTouchInputModuleEditor.cs ver. 1.0 */
using UnityEngine;
using UnityEditor;
using System.Collections;
#pragma warning disable 618

namespace UnityEngine.EventSystems
{
	[CustomEditor( typeof( StandaloneMultiTouchInputModule ) )]
	public class StandaloneMultiTouchInputModuleEditor : Editor
	{
		public void OnEnable ()
		{
			// Store the joystick that we are selecting
			StandaloneMultiTouchInputModule smtim = ( StandaloneMultiTouchInputModule )target;

			if( smtim.gameObject.GetComponent<StandaloneInputModule>() && smtim.gameObject.GetComponent<StandaloneInputModule>().enabled == true )
			{
				Debug.Log( "'Standalone Input Module' has been disabled so that the 'Standalone Multi-Touch Input Module' can work correctly." );
				smtim.gameObject.GetComponent<StandaloneInputModule>().enabled = false;
			}
			if( smtim.gameObject.GetComponent<TouchInputModule>() && smtim.gameObject.GetComponent<TouchInputModule>().enabled == true )
			{
				Debug.Log( "'Touch Input Module' has been disabled so that the 'Standalone Multi-Touch Input Module' can work correctly." );
				smtim.gameObject.GetComponent<TouchInputModule>().enabled = false;
			}
		}

		public override void OnInspectorGUI ()
		{
			// Store the joystick that we are selecting
			StandaloneMultiTouchInputModule smtim = ( StandaloneMultiTouchInputModule )target;

			if( smtim.gameObject.GetComponent<StandaloneInputModule>() && smtim.gameObject.GetComponent<StandaloneInputModule>().enabled == true )
			{
				Debug.Log( "'Standalone Input Module' has been disabled so that the 'Standalone Multi-Touch Input Module' can work correctly." );
				smtim.gameObject.GetComponent<StandaloneInputModule>().enabled = false;
			}
			if( smtim.gameObject.GetComponent<TouchInputModule>() && smtim.gameObject.GetComponent<TouchInputModule>().enabled == true )
			{
				Debug.Log( "'Touch Input Module' has been disabled so that the 'Standalone Multi-Touch Input Module' can work correctly." );
				smtim.gameObject.GetComponent<TouchInputModule>().enabled = false;
			}
		}
	}
}