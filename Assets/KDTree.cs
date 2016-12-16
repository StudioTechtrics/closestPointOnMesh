//#define DEBUG_KD
using UnityEngine;
using System.Collections;
using System.Collections.Generic;
/// <summary>
/// Recursively paritions a mesh's vertices to allow to more quickly
/// narrow down the search for a nearest point on it's surface with respect to another
/// point
/// </summary>
/// 
using UniRx;

[RequireComponent(typeof(MeshCollider))]
public class KDTree : MonoBehaviour {

    static Dictionary<int, KDTreeData> dict = new Dictionary<int, KDTreeData>();

    public System.Action OnComputationFinish = null;

    KDTreeData data = null;

    public bool autoInit = false;

    void Awake() { if (autoInit) Init(); }

    public bool Processing {
        get {

            if (data == null) return true;
            return data.Processing;
        }
        private set {

            if (data == null) return;
            data.Processing = value;
        }
    }

    void Init() {

        var sharedMesh = GetComponent<MeshCollider>().sharedMesh;

        // Already contains
        if (dict.ContainsKey(sharedMesh.GetInstanceID())) {

            // if already contains, do nothing!!!
            data = dict[sharedMesh.GetInstanceID()];

            if (OnComputationFinish != null)
                OnComputationFinish();
        }
        else {

            data = new KDTreeData(sharedMesh);

            //Debug.Log("Generating tree for " + data.ID);

            dict[data.ID] = data;

            // build kd-tree
            Processing = true;
            var backgroundMethod = Observable.Start(() => {

                data.Build();
            });

            Observable
            .WhenAll(backgroundMethod)
            .ObserveOnMainThread()
            .Subscribe((unit) => {

                Processing = false;
                if (OnComputationFinish != null)
                    OnComputationFinish();       
            });
        }
    }

    public bool ClosestPointOnOptimized(Vector3 to, out Vector3 closestPoint, /*ref int triangleCacheIndex,*/ float range = -1f) {
        if (Processing)
        {   //Check yourself before you wreck yourself
            Debug.LogError("Wait before KD-Tree has been built before using it");

            closestPoint = to;
            return false;
        }

        return data.ClosestPointOnOptimized(to, out closestPoint, transform, range);
    }

    public bool visualization = false;

    void OnDrawGizmos() {

        if (!visualization)
            return;

        Gizmos.matrix = transform.localToWorldMatrix;
        
        if(data != null && !Processing)
            data.RecursiveDraw(data.rootNode, 0, true);
    }
}