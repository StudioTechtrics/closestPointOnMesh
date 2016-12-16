using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(LineRenderer))]
public class ClosestPointOnMeshTest : MonoBehaviour {

    [Range(0.01f, 100f)]
    public float QueryRadius = 15f;

	void Start () {
	    	
	}
	
	void Update () {

        Vector3 pos = transform.position;

        Collider[] colls = Physics.OverlapSphere(pos, QueryRadius);


	}
}