"""Library of handy plotting functions."""
import numpy as np
import matplotlib.pyplot as plt


def get_MQ_UQ_d(bam_df):
  mq = np.concatenate([bam_df['m1_MQ'].values, bam_df['m2_MQ'].values])
  idx = (mq < 255)
  mq = mq[idx]
  uq = np.concatenate([bam_df['m1_UQ'].values, bam_df['m2_UQ'].values])[idx]
  d  = np.concatenate([bam_df['m1_d_error'].values, bam_df['m2_d_error'].values])[idx]
  return mq, uq, d


def plot_MQ_vs_wrong(bam_df, t=10, fmt='ko', label='line'):
  """Plots fraction of incorrect reads (as defined by d-threshold value) on a log axis, against MQ value.
  The theoretical curve for this quantity is also plotted."""
  mq, uq, d = get_MQ_UQ_d(bam_df)
  idx_1 = (d < t)
  idx_2 = (d >= t)
  mq_bins = np.arange(-0.5, 70.5, step=1.0)
  mq_bin_centers = (mq_bins[:1] + mq_bins[1:])/2.0
  # uq_bins = np.arange(-0.5, 100.5, step=1.0)

  mq_correct, _ = np.histogram(mq[idx_1], bins=mq_bins)
  mq_wrong, _ = np.histogram(mq[idx_2], bins=mq_bins)

  #f = np.clip(mq_wrong.astype(float) / (mq_correct + mq_wrong), a_min=1e-6, a_max=2)
  f = mq_wrong.astype(float) / (mq_correct + mq_wrong)
  f_theoretical = 10.0**(-mq_bin_centers/10.0)

  plt.semilogy(mq_bin_centers, f, fmt, label=label)
  plt.semilogy(mq_bin_centers, f_theoretical, 'k:')
  plt.xlabel('MQ')
  plt.ylabel(r'$p_{wrong}$')


def plot_MQ_balance(bam_df, t=10, fmt='k', label='line'):
  """Given a d-threshold for indicating correct/incorrect give us percentage
  of total number of wrong reads greater than given MQ"""
  mq, uq, d = get_MQ_UQ_d(bam_df)
  idx_2 = (d >= t)
  mq_bins = np.arange(-0.5, 70.5, step=1.0)
  mq_bin_centers = (mq_bins[:1] + mq_bins[1:])/2.0

  mq_wrong, _ = np.histogram(mq[idx_2], bins=mq_bins)

  f = mq_wrong[::-1].cumsum()[::-1]
  plt.step(mq_bin_centers, f / float(mq.size), fmt, label=label, drawstyle='steps-post')
  plt.gca().set_yscale('log')
  plt.xlabel('MQ')
  plt.ylabel('Fraction\nreads wrong')
  plt.title('Fraction of reads wrong for MQ >= value')


def plot_UQ_vs_wrong(bam_df, t=10, fmt='ko', label='line'):
  """Plots fraction of incorrect reads (as defined by d-threshold value) on a log axis, against UQ value."""
  mq, uq, d = get_MQ_UQ_d(bam_df)
  idx_1 = (d < t)
  idx_2 = (d >= t)
  uq_bins = np.arange(-0.5, 100.5, step=1.0)
  uq_bin_centers = (uq_bins[:1] + uq_bins[1:])/2.0

  uq_correct, _ = np.histogram(uq[idx_1], bins=uq_bins)
  uq_wrong, _ = np.histogram(uq[idx_2], bins=uq_bins)

  f = uq_wrong.astype(float) / (uq_correct + uq_wrong)

  plt.semilogy(uq_bin_centers, f, fmt, label=label)
  plt.xlabel('UQ')
  plt.ylabel(r'$p_{wrong}$')


def plot_UQ_balance(bam_df, t=10, fmt='k', label='line'):
  """Given a d-threshold for indicating correct/incorrect give us percentage
  of total number of wrong reads greater than given UQ"""
  mq, uq, d = get_MQ_UQ_d(bam_df)
  idx_2 = (d >= t)
  uq_bins = np.arange(-0.5, 100.5, step=1.0)
  uq_bin_centers = (uq_bins[:1] + uq_bins[1:])/2.0

  mq_wrong, _ = np.histogram(mq[idx_2], bins=uq_bins)

  f = mq_wrong[::-1].cumsum()[::-1]
  plt.step(uq_bin_centers, f / float(mq.size), fmt, label=label, drawstyle='steps-post')
  plt.gca().set_yscale('log')
  plt.xlabel('UQ')
  plt.ylabel('Fraction\nreads wrong')
  plt.title('Fraction of reads wrong for UQ >= value')




def plot_MQ_vs_UQ(bam_df, t = 10):
  mq = np.concatenate([bam_df['m1_MQ'].values, bam_df['m2_MQ'].values])
  idx = (mq < 255)
  mq = mq[idx]
  uq = np.concatenate([bam_df['m1_UQ'].values, bam_df['m2_UQ'].values])[idx]
  d  = np.concatenate([bam_df['m1_d_error'].values, bam_df['m2_d_error'].values])[idx]
  idx_1 = (d < t)
  idx_2 = (d >= t)
  plt.subplot(1, 2, 1)
  plt.plot(mq[idx_1], uq[idx_1], 'g.')
  plt.subplot(1, 2, 2)
  plt.plot(mq[idx_2], uq[idx_2], 'r.')
  plt.xlabel('MQ')
  plt.ylabel('UQ')
  plt.title('MQ vs UQ at t = {}'.format(t))


def plot_MQ_UQ_hist(bam_df, t = 10):
  """."""
  mq = np.concatenate([bam_df['m1_MQ'].values, bam_df['m2_MQ'].values])
  uq = np.concatenate([bam_df['m1_UQ'].values, bam_df['m2_UQ'].values])
  d  = np.concatenate([bam_df['m1_d_error'].values, bam_df['m2_d_error'].values])
  idx_1 = (d < t)
  idx_2 = (d >= t)
  plt.subplot(2, 2, 1)
  plt.hist(mq[idx_1], uq[idx_1], 'g.')
  plt.subplot(1, 2, 2)
  plt.plot(mq[idx_2], uq[idx_2], 'r.')
  plt.xlabel('MQ')
  plt.ylabel('UQ')
  plt.title('MQ vs UQ at t = {}'.format(t))
