Shader "Unlit/Line"
{
    Properties
    {
        _LineWight ("Line Wight", Float) = 0.1
    }
    SubShader
    {
        Tags { "Queue"="Transparent" "IgnoreProjector"="True" "RenderType"="Transparent" "PreviewType"="Plane" "DisableBatching"="True"}
        Cull Off
        Blend SrcAlpha OneMinusSrcAlpha
        ZWrite Off

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            // make fog work
            #pragma multi_compile_fog

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float4 color : COLOR;
                float3 normal : NORMAL;
                float2 texcoord : TEXCOORD0;
            };

            struct v2f
            {
                // float2 uv : TEXCOORD0;
                float4 vertex : SV_POSITION;
                float4 color : COLOR;
                // float3 normal : NORMAL;
            };


            float _LineWight;

            v2f vert (appdata v)
            {
                v2f o;
                o.color = v.color;
                // o.normal = v.normal;
                // o.uv = v.texcoord;

                float3 dir2cam = ObjSpaceViewDir(v.vertex);
                dir2cam = normalize(dir2cam);

                v.normal = mul((float3x3)unity_WorldToObject, v.normal);
                float3 pushDir = cross(v.normal, dir2cam)*v.texcoord.x;
                pushDir = normalize(pushDir);

                v.vertex.xyz+=pushDir*_LineWight;

                o.vertex = UnityObjectToClipPos(v.vertex);

                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                return i.color;
                // return float4(i.normal, 1);
            }
            ENDCG
        }
    }
}
