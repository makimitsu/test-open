import tensorflow as tf
import scipy.io
import os
import pathlib
import sys
import numpy as np

date = '240621' #sys.argv[1]
shot = '41' #sys.argv[2]
N_projection = 50 # int(sys.argv[3])
N_grid = 91 #int(sys.argv[4])

PATH = pathlib.Path('/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data/SXR_data/' + date + '/shot' + shot)
savePATH = '/Users/shohgookazaki/Library/CloudStorage/GoogleDrive-shohgo-okazaki@g.ecc.u-tokyo.ac.jp/My Drive/OnoLab/data/result_matrix'

def get_sorted_file_list(dataset_path, pattern):
    file_list = list(dataset_path.glob(pattern))
    file_list.sort()  # ファイル名順にソート
    return file_list

def load_mat_files(input_path):
    mat_sxr = scipy.io.loadmat(input_path.numpy())
    return mat_sxr['sxr1'], mat_sxr['sxr2'],mat_sxr['sxr3'],mat_sxr['sxr4']  # 必要に応じてキーを変更してください

def load_data(input_data):
    sxr1, sxr2, sxr3, sxr4 = load_mat_files(input_data)
    
    sxr1 = tf.expand_dims(sxr1, -1)
    sxr2 = tf.expand_dims(sxr2, -1)
    sxr3 = tf.expand_dims(sxr3, -1)
    sxr4 = tf.expand_dims(sxr4, -1)
    
    sxr1 = tf.cast(sxr1, tf.float32)
    sxr2 = tf.cast(sxr2, tf.float32)
    sxr3 = tf.cast(sxr3, tf.float32)
    sxr4 = tf.cast(sxr4, tf.float32)
    
    sxr1 = tf.image.flip_up_down(sxr1)
    sxr2 = tf.image.flip_up_down(sxr2)
    sxr3 = tf.image.flip_up_down(sxr3)
    sxr4 = tf.image.flip_up_down(sxr4)

    return sxr1, sxr2, sxr3, sxr4


def tf_load_data(input_data):
    return tf.py_function(func=load_data, inp=[input_data], Tout=[tf.float32, tf.float32, tf.float32, tf.float32])

def downsample(filters, size, apply_batchnorm=True):
  initializer = tf.random_normal_initializer(0., 0.02)

  result = tf.keras.Sequential()
  result.add(
      tf.keras.layers.Conv2D(filters, size, strides=2, padding='same',
                             kernel_initializer=initializer, use_bias=False))

  if apply_batchnorm:
    result.add(tf.keras.layers.BatchNormalization())

  result.add(tf.keras.layers.LeakyReLU())

  return result

def upsample(filters, size, apply_dropout=False):
  initializer = tf.random_normal_initializer(0., 0.02)

  result = tf.keras.Sequential()
  result.add(
    tf.keras.layers.Conv2DTranspose(filters, size, strides=2,
                                    padding='same',
                                    kernel_initializer=initializer,
                                    use_bias=False))

  result.add(tf.keras.layers.BatchNormalization())

  if apply_dropout:
      result.add(tf.keras.layers.Dropout(0.5))

  result.add(tf.keras.layers.ReLU())

  return result

def Generator():
  inputs = tf.keras.layers.Input(shape=[N_projection,N_projection,1])

  down_stack = [
    downsample(64, 4, apply_batchnorm=False),  # (batch_size, 128, 128, 64)
    downsample(128, 4),  # (batch_size, 64, 64, 128)
    downsample(256, 4),  # (batch_size, 32, 32, 256)
    downsample(512, 4),  # (batch_size, 16, 16, 512)
    downsample(512, 4),  # (batch_size, 8, 8, 512)
    downsample(512, 4),  # (batch_size, 4, 4, 512)
    #downsample(512, 4),  # (batch_size, 2, 2, 512)
    #downsample(512, 4),  # (batch_size, 1, 1, 512)
  ]

  up_stack = [
    #upsample(512, 4, apply_dropout=True),  # (batch_size, 2, 2, 1024)
    upsample(512, 4, apply_dropout=True),  # (batch_size, 4, 4, 1024)
    upsample(512, 4, apply_dropout=True),  # (batch_size, 8, 8, 1024)
    upsample(512, 4),  # (batch_size, 16, 16, 1024)
    upsample(256, 4),  # (batch_size, 32, 32, 512)
    upsample(128, 4),  # (batch_size, 64, 64, 256)
    upsample(64, 4),  # (batch_size, 128, 128, 128)
  ]

  initializer = tf.random_normal_initializer(0., 0.02)
  
  last = tf.keras.layers.Conv2DTranspose(OUTPUT_CHANNELS, 4,
                                         strides=1,
                                         padding='same',
                                         kernel_initializer=initializer,
                                         activation='tanh')  # (batch_size, 256, 256, 3)
  
  x = inputs

  # Downsampling through the model
  skips = []
  for down in down_stack:
    x = down(x)
    skips.append(x)

  skips = reversed(skips[:-1])
  
  
  for up, skip in zip(up_stack, skips):
    x = up(x)
    x = tf.keras.layers.Lambda(lambda inputs: tf.image.resize(inputs[0], size=(tf.shape(inputs[1])[1], tf.shape(inputs[1])[2])))([x, skip])
    x = tf.keras.layers.Concatenate()([x, skip])
  
  x = tf.keras.layers.Lambda(lambda x: tf.image.resize(x, size=(N_grid, N_grid)))(x)
  x = last(x)

  return tf.keras.Model(inputs=inputs, outputs=x)

def Discriminator():
  initializer = tf.random_normal_initializer(0., 0.02)

  inp = tf.keras.layers.Input(shape=[N_grid,N_grid,1], name='input_image')
  tar = tf.keras.layers.Input(shape=[N_grid,N_grid,1], name='target_image')

  x = tf.keras.layers.concatenate([inp, tar])  # (batch_size, 256, 256, channels*2)

  down1 = downsample(64, 4, False)(x)  # (batch_size, 128, 128, 64)
  down2 = downsample(128, 4)(down1)  # (batch_size, 64, 64, 128)
  down3 = downsample(256, 4)(down2)  # (batch_size, 32, 32, 256)

  zero_pad1 = tf.keras.layers.ZeroPadding2D()(down2)  # (batch_size, 34, 34, 256)
  conv = tf.keras.layers.Conv2D(512, 4, strides=1,
                                kernel_initializer=initializer,
                                use_bias=False)(zero_pad1)  # (batch_size, 31, 31, 512)

  batchnorm1 = tf.keras.layers.BatchNormalization()(conv)

  leaky_relu = tf.keras.layers.LeakyReLU()(batchnorm1)

  zero_pad2 = tf.keras.layers.ZeroPadding2D()(leaky_relu)  # (batch_size, 33, 33, 512)

  last = tf.keras.layers.Conv2D(1, 4, strides=1,
                                kernel_initializer=initializer)(zero_pad2)  # (batch_size, 30, 30, 1)

  return tf.keras.Model(inputs=[inp, tar], outputs=last)

def save_images(model, index, sxr1, sxr2, sxr3, sxr4):
  
  EE1 = np.squeeze(model(sxr1, training=True))
  EE2 = np.squeeze(model(sxr2, training=True))
  EE3 = np.squeeze(model(sxr3, training=True))
  EE4 = np.squeeze(model(sxr4, training=True))
  
  EE1 = tf.cast(EE1, tf.float64)
  EE2 = tf.cast(EE2, tf.float64)
  EE3 = tf.cast(EE3, tf.float64)
  EE4 = tf.cast(EE4, tf.float64)
     
  
  # Create a dictionary to hold the data
  data_dict = {
      'EE1': EE1,
      'EE2': EE2,
      'EE3': EE3,
      'EE4': EE4
  }

  # Save the data to a .mat file
  savedir = savePATH + '/cGAN/' + date + '/' + 'shot' + shot 
  
  if not os.path.exists(savedir):
    os.makedirs(savedir)
  savefile = savedir + '/' + f'{index}.mat'
  scipy.io.savemat(savefile, data_dict)


# The facade training set consist of 400 images
BUFFER_SIZE = 800
# The batch size of 1 produced better results for the U-Net in the original pix2pix experiment
BATCH_SIZE = 1

input_files = get_sorted_file_list(PATH, '*.mat')

input_dataset = tf.data.Dataset.from_tensor_slices([str(f) for f in input_files])
input_dataset = input_dataset.map(tf_load_data, num_parallel_calls=tf.data.AUTOTUNE)
input_dataset = input_dataset.batch(BATCH_SIZE)

OUTPUT_CHANNELS = 1

generator = Generator()


discriminator = Discriminator()

generator_optimizer = tf.keras.optimizers.Adam(2e-5, beta_1=0.5)
discriminator_optimizer = tf.keras.optimizers.Adam(2e-5, beta_1=0.5)

checkpoint_dir = './training_checkpoints'
checkpoint_prefix = os.path.join(checkpoint_dir, "ckpt")
checkpoint = tf.train.Checkpoint(generator_optimizer=generator_optimizer,
                                 discriminator_optimizer=discriminator_optimizer,
                                 generator=generator,
                                 discriminator=discriminator)

  
# Restoring the latest checkpoint in checkpoint_dir
checkpoint.restore(tf.train.latest_checkpoint(checkpoint_dir))

# Run the trained model on a few examples from the test set
for index, (sxr1, sxr2, sxr3, sxr4) in enumerate(input_dataset):
  save_images(generator,index+1, sxr1, sxr2, sxr3, sxr4)