using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;


namespace UnityStandardAssets.Vehicles.Car
{
    [RequireComponent(typeof(CarController))]
    public class SphereAI4 : MonoBehaviour
    {

        List<GameObject> FindGameObjectWithLayer(int layer)
        {
            List<GameObject> good_list = new List<GameObject>();
            var allObjects = FindObjectsOfType<GameObject>();
            foreach(GameObject objectToCheck in allObjects){
                if(objectToCheck.layer == layer)
                {
                    good_list.Add(objectToCheck);
                }
            }
            return good_list;
        }

        private CarController m_Car; // the car controller we want to use

        public GameObject terrain_manager_game_object;
        TerrainManager terrain_manager;

        GameObject leader;
        const int num_followers = 4;
        Vector3[] targets = new Vector3[num_followers];
        Vector3[] target_vels = new Vector3[num_followers];
        Vector3[] friends_pos = new Vector3[num_followers];
        Vector3[] friend_vels = new Vector3[num_followers];
        int[] offsets = new int[]{-1, 1, -2, 2};
        Vector3 leader_pos;
        Vector3 leader_back;
        Vector3 leader_left;
        int my_target;
        Vector3 my_old_position;
        bool backup;
        GameObject[] friends;

        GameObject mySphere;
        Vector3 magic_pos;
        int counter = 0;

        private void Start()
        {
            friends = GameObject.FindGameObjectsWithTag("Player");
            // get the car controller
            m_Car = GetComponent<CarController>();
            terrain_manager = terrain_manager_game_object.GetComponent<TerrainManager>();

            leader = GameObject.Find("ReplayCar (2)");
            
            switch(transform.position.x) {
                case 239f:
                    my_target = 0;
                    break;
                case 258.3f:
                    my_target = 1;
                    break;
                case 246.8f:
                    my_target = 2;
                    break;
                case 252.6f:
                    my_target = 3;
                    break;
                default:
                    Debug.Log("didn't find target");
                    break;
            }

            my_old_position = transform.position;
        }

        //(2) = 239, 230
        //(3) = 258.3, 230
        //(4) = 246.8, 230
        //(5) = 252.6, 230


        private void FixedUpdate()
        {
            if (counter < 40)
            {
                counter++;
                for (int i = 0; i < num_followers; i++)
                {
                    targets[i] = leader.transform.position;
                }
                return;
            }

            leader_pos = leader.transform.position;
            leader_back = leader.transform.TransformDirection(Vector3.back);
            leader_left = leader.transform.TransformDirection(Vector3.left);
            for (int i = 0; i < num_followers; i++) {
                Vector3 new_pos = leader_pos /*+ leader_back * 10f * (i + 1) + leader_left * offsets[i] * 20f*/;
                Vector3 diff = (new_pos - targets[i]);
                target_vels[i] = (new_pos - targets[i]) / Time.fixedDeltaTime;
                targets[i] = new_pos;
            }
            foreach(var target in targets) {
                Debug.DrawLine(leader.transform.position, target, Color.red);
            }

            Vector3 q_tar = targets[my_target];
            Vector3 v_tar = target_vels[my_target];
            float theta_tar = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, v_tar, Vector3.down);
            Vector3 q_obs = Vector3.positiveInfinity;

           //foreach (GameObject obstacle in GameObject.FindGameObjectsWithTag("cube"))
           foreach (GameObject obstacle in FindGameObjectWithLayer(9))
           {
                Vector3 closest_point = obstacle.GetComponent<Collider>().ClosestPoint(transform.position);
                Vector3 min_distance = closest_point - transform.position;    
                if (min_distance.sqrMagnitude < (q_obs - transform.position).sqrMagnitude)
                {
                    q_obs = closest_point;
                }
            }

            Vector3 q_ao = q_obs - transform.position;

            /*for(int i = 0; i < num_followers; i++)
            {
                Vector3 new_friend_pos = friends[i].transform.position;
                friend_vels[i] = (new_friend_pos - friends_pos[i]) * Time.fixedDeltaTime;
                friends_pos[i] = new_friend_pos;
            }*/

            Vector3 q_at = q_tar - transform.position;
            float psi = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_at, Vector3.down);
            float xi_1 = 2f;
            float xi_2 = 1f;
            float rho_inv = 1f / (q_obs - transform.position).magnitude;
            float rho_zero_inv = 1f / 20f;
            float mu = (xi_2 * rho_inv * (rho_inv - rho_zero_inv)) / q_ao.magnitude;
            float lambda = (mu * q_ao.magnitude) / (xi_1 * q_at.magnitude);
            float theta_ao = Mathf.Deg2Rad * Vector3.SignedAngle(Vector3.right, q_ao, Vector3.down);
            float psi_bar = Mathf.Atan((Mathf.Sin(psi) - lambda * Mathf.Sin(theta_ao)) / (Mathf.Cos(psi) - lambda * Mathf.Cos(theta_ao)));
            float term_1 = v_tar.magnitude * Mathf.Cos(theta_tar - psi);
            //Debug.Log("Term 1: " + term_1);
            float term_2 = xi_1 * q_at.magnitude;
            //Debug.Log("Term 2: " + term_2);
            float term_3 = v_tar.sqrMagnitude * Mathf.Pow(Mathf.Sin(theta_tar - psi_bar), 2);
            //Debug.Log("Term 3: " + term_3);
            float norm_v_i = Mathf.Sqrt(Mathf.Pow(v_tar.magnitude * Mathf.Cos(theta_tar - psi) + xi_1 * q_at.magnitude, 2) + v_tar.sqrMagnitude * Mathf.Pow(Mathf.Sin(theta_tar - psi_bar), 2));
            float theta_i = theta_tar;
            if (norm_v_i != 0)
            {
                theta_i = psi_bar + Mathf.Asin(Mathf.Clamp((v_tar.magnitude * Mathf.Sin(theta_tar - psi_bar) / norm_v_i), -1, 1));
            }
            Debug.Log("psi_bar " + psi_bar * Mathf.Rad2Deg);
            Debug.Log("theta_ao " + theta_ao);
            Debug.Log("psi " + psi);

            Vector3 u_i = new Vector3(norm_v_i * -Mathf.Cos(theta_i), 0f, norm_v_i * -Mathf.Sin(theta_i));

            my_old_position = transform.position;
            transform.position += u_i * Time.fixedDeltaTime;

            //Debug.DrawLine(targets[my_target], targets[my_target] + new Vector3(Mathf.Cos(theta_tar), 0, Mathf.Sin(theta_tar)) * 4f);
            Debug.DrawLine(targets[my_target], targets[my_target] + v_tar, Color.blue);
            Debug.DrawLine(transform.position, transform.position + new Vector3(Mathf.Cos(theta_i), 0, Mathf.Sin(theta_i)) * 10f);
            Debug.DrawLine(transform.position, transform.position + new Vector3(Mathf.Cos(theta_ao), 0, Mathf.Sin(theta_ao)) * 10f, Color.green);
            Debug.DrawLine(transform.position, transform.position + new Vector3(Mathf.Cos(psi), 0, Mathf.Sin(psi)) * 10f, Color.cyan);
        }
    }
}
