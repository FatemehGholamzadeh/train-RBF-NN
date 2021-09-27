# train-RBF-NN
training a RBF network Using ES algorithm - Computational Intelligence Course (Spring 2019)  

in this project we implement Evolution Strategy Algorithm and use it to train a RBF Neural Network.  
  
  


<img width="425" alt="rbf" src="https://user-images.githubusercontent.com/44861408/134909024-5cdede84-5553-4703-8965-bc01a5830d47.png">  


we store outputs of first layer in G Matrix: 

![image](https://user-images.githubusercontent.com/44861408/134909701-f43a66e2-4829-4183-ad6c-2c16b656792c.png)
  
    
    
by performing linear transformation output can be written as : 

![image](https://user-images.githubusercontent.com/44861408/134909834-37317218-7db1-4398-b443-94ee04dbc6d7.png)

![image](https://user-images.githubusercontent.com/44861408/134909877-4e2a4f0d-d35f-4172-b5bc-c59cf860f95d.png)  
  
  
 
so the Loss Function can be define as bellow: 

![image](https://user-images.githubusercontent.com/44861408/134908552-76766568-df13-4a11-893d-46ff958c8454.png)  
  
  
and weights can be calculated like this : 

![image](https://user-images.githubusercontent.com/44861408/134910009-76a2eb01-13c9-42f6-857f-e070f88529e1.png)  
  
  

then we use this network for classification and regression tasks.

## Regression Task

the output for Regression Task (on regdata2000.xhxs) is show in bellow:

![image](https://user-images.githubusercontent.com/44861408/134909361-29d11c30-ee1b-4a07-a41f-612ff87dfc12.png)

